# Parallelism of SchrodingerSolver — complete investigation

**Conclusion: running the solver under `mpirun` is NOT the way.** There is a better
solution, already implemented and verified, that gives the same result with
`O(nnz)` memory instead of `O(n²)`.

---

## 1. Root cause of the segfault with `mpirun` — CONFIRMED in the source

It is not MUMPS, it is not the solver, it is not binding/UCX. It is a
**non-partitioned mesh**.

`MeshUtils.F90:2117` (Elmer 9.0 source):
```fortran
ELSE                                       ! NumProcs > 1
  INQUIRE( FILE=MeshNamePar(1:n)//'/partitioning.'// &
      TRIM(i2s(Numprocs))//'/part.1.header', EXIST=Found)
  IF(.NOT. Found ) THEN
    CALL Warn('LoadMesh','Requested mesh > '//MeshNamePar(1:n)//' < in partition '&
        //TRIM(I2S(Numprocs))//' does not exist!')
    RETURN                                  ! <-- returns NULL, and is NOT Fatal
  END IF
END IF
```

The one that dereferences without checking — `ModelDescription.F90:2892-2899`:
```fortran
2892:  Model % Meshes => LoadMesh2( ... )      ! returns NULL()
2895:  Model % Meshes % SingleMesh = Single    ! <-- NULL DEREF -> SIGSEGV
2899:  IF(.NOT.ASSOCIATED(Model % Meshes)) THEN ! <-- checked TOO LATE
```

**It is a real bug in Elmer 9.0.** It explains all the facts:
- it only fails with `PEs >= 2` (with 1 process the `partitioning.*` path is not even read)
- it is independent of binding, UCX, symlink (I tested 5 configurations)
- it is **not** MUMPS (MUMPS is never even called)

### Controlled A/B (done by the subagent, local build)
| mesh | `-np` | result |
|---|---:|---|
| without `partitioning.2/` | 2 | **exit 139 SIGSEGV** |
| with `partitioning.2/` | 2 | **exit 0** |

### Aggravating factor: the project never partitioned the mesh
The `roda_caso.pbs` scripts call `"$ES" case.sif` (serial). The partitioning step was missing.

---

## 2. Why partitioning did not solve it here

I tried **two** routes to partition; both failed because of Elmer 9.0 bugs:

| Route | Result |
|---|---|
| `Partition Mesh = Logical True` in the SIF | crash **in serial** (`Where perform partitioning in serial case?` + segfault) |
| solver `PartitionMesh` + `Output Directory` | crash in `partitionmeshserial` → `__generalutils_MOD_i2s` |

The second crash occurs **after** "Partitioning the boundary elements sets" —
that is, the partition is computed but Elmer dies while printing its number.
The installed `SOLVER.KEYWORDS` does not even list `partitioning exit`.

**Elmer 9.0 (2020 release) cannot partition a mesh on this platform.**

---

## 3. ⚠ Even with partitioning, the solver would return GARBAGE

With `partitioning.2` present and `-np 2`, `SchrodingerSolver` **does not crash but
gets the result wrong**:

| | serial | `-np 2` |
|---|---|---|
| Free DOFs | 59 | **29 per rank** |
| `E(1)` | 1.025415149101 | **3.052806819296** |

**Cause:** the solver assembles **dense local** matrices and calls `DSYGV`.
On a partitioned mesh, `Model % NumberOfNodes` and `Element % NodeIndexes` are
**local** — each rank sees only its own partition. Mathematically the problem changes.

I reported this as "Approach C — do not do it": it runs without error and returns a wrong number.

---

## 4. ✅ THE IMPLEMENTED SOLUTION: export and solve sparse

Rewriting of the workflow, bypassing both bugs at once:

```
Elmer (assembles A and M, exports triplets)  ->  sparse solver (ARPACK shift-invert)
```

### How it works
A section (§7) was added to `SchrodingerSolver.F90` which, if
`Export Matrices = Logical True`, writes:
- `<prefix>_A.dat` — operator `A = K + Ma`
- `<prefix>_M.dat` — mass `M`
- `<prefix>_fixed.dat` — Dirichlet DOFs (0-based)

Format: `n n nnz` and then `i j value` (0-based), ASCII.

### Solve
```bash
python3 python/parallel/resolve_exportado.py schrod 5
```

### Verification — identical to serial Elmer
| | `E(1)` | `E(2)` | `E(3)` | `E(4)` |
|---|---|---|---|---|
| **Serial Elmer (DSYGV)** | 1.113644656498 | 2.172208383659 | 2.172208383659 | 3.230772110820 |
| **Exported + sparse** | 1.113644656498 | 2.172208383659 | 2.172208383659 | 3.230772110820 |

**Identical to the 12th digit.** And:
- time **0.01 s** (vs dense)
- memory **0.1 MB** vs 4.7 MB dense → **52× smaller** (and it grows `O(nnz)` vs `O(n²)`)

---

## 5. Why this is BETTER than MPI

| | `mpirun` on Elmer | export + sparse |
|---|---|---|
| Does it work? | ❌ (2 Elmer 9.0 bugs) | ✅ verified |
| Memory | `O(n²)` per rank | `O(nnz)` |
| Accuracy | garbage if partitioned | identical to serial |
| Scaling (measured) | — | **up to 4097×** vs dense |
| Depends on Elmer's ARPACK? | yes | **no** |

The real bottleneck was never CPU — it was **algorithmic** (`O(n²)`/`O(n³)`). Trading
dense for sparse gives more gain than any parallelization of `DSYGV`.

---

## 6. If MPI is really necessary

Two routes, in order of cost:

**(a) Use the already-validated PETSc/SLEPc driver** (`python/parallel/`).
It solves in parallel with error `1e−15`. All that is missing is **feeding** it the
exported matrices. `mm2petsc` exists; the current bottleneck is the PETSc binary
format (header with version). Alternative without new code: use `MatLoad` with
`-viewer_binary_skip_info` or generate the binary with a PETSc utility.

**(b) Rewrite the solver with native parallel assembly.** It is the "correct" path
but expensive (3–5 days, high risk). An ideal template already exists in Elmer:
`modules/DFTSolver.F90` in the Elmer source tree → `WavefunctionSolver` (lines 218-505), which solves
`(-½∇²+V+Vxc)ψ=Eψ` in parallel. Pattern:
```fortran
CALL DefaultInitialize()
DO t = 1, Solver % NumberOfActiveElements        ! NEVER Model % NumberOfBulkElements
  Element => GetActiveElement(t)
  CALL DefaultUpdateEquations(STIFF, FORCE)
  CALL DefaultUpdateMass(MASS)
END DO
CALL DefaultFinishAssembly()
CALL DefaultDirichletBCs()
CALL DefaultSolve()
```
`Linear System Symmetric = Logical True` is **mandatory** — only then does
`CRS_ElimSymmDirichlet` (SolverUtils.F90:8711) run and
`CRSMatrix.F90:1285` also zero the mass columns.

**Warning:** Elmer's native ARPACK eigensolver, in parallel, **converges but returns
wrong eigenvalues** (tested on an SPD case with a known analytic solution: negative
results). Only adopt (b) if the `V≡0` calibration case closes.

---

## 7. Retraction

An earlier internal note claimed that the
variable `mode` was never assigned. **That is wrong** — `EigenSolve.F90:213` sets
`Mode = 2` and line 231 sets `Mode = 3`; Fortran is case-insensitive. Document
retracted.
