# schrodinger_P1_2 — Schrödinger Solver in Elmer FEM

**Status: WORKING AND VALIDATED in 1D, 2D and 3D.**
Elmer **does not have** a Schrödinger solver; this one was written from scratch.

---

## 1. Why a solver of its own

Elmer 9.0 has two defects in its ARPACK eigensolver:

1. **`EigenSolve.F90` does not initialize `RESID`** (the ARPACK starting vector) in
   `ArpackEigenSolve` → `info = -9` (*starting vector is zero*).
   Root cause diagnosed and documented in `../schrodinger_P1/RELATORIO_OPCAO_B.md`.
2. Even with `RESID` fixed (patch applied), the *reverse communication* loop
   fails with `-9` after the first application of the operator — a second, distinct defect.

**Solution:** a dedicated solver that assembles the matrices and calls **`DSYGV` (LAPACK)**
directly, without ARPACK.

---

## 2. Formulation

Equation (ħ = m = 1): `-∇²ψ + V(x)ψ = Eψ`

Weak Galerkin form → generalized eigenproblem `A ψ = E M ψ`, `A = K + Ma`:

| Matrix | Integral | Assembly in the code |
|---|---|---|
| `K` stiffness | `∫∇Nᵢ·∇Nⱼ dΩ` | `Σ s·detJ·(dBasisdx(p,:)·dBasisdx(q,:))` |
| `M` mass | `∫NᵢNⱼ dΩ` | `Σ s·detJ·Basis(p)·Basis(q)` |
| `Ma` potential | `∫Nᵢ V Nⱼ dΩ` | `Σ s·detJ·V(ξ)·Basis(p)·Basis(q)` |

Quadrature via `GaussPoints(Element)`; `detJ` and the shape functions via `ElementInfo`.
Generic dimensionality (1D/2D/3D) via `CoordinateSystemDimension()` and
`dBasisdx(p,1:dim)`.

---

## 3. Validation — cross-verification in three dimensions

Exact analytic solution in all cases. Python and Elmer have **independent assembly**
and a **different eigensolver** (scipy `eigh` vs LAPACK `DSYGV`).

### 1D — infinite well, L=40, 250 elements `Eₙ=(nπ/L)²`

| n | Elmer | analytic | error |
|---|---|---|---|
| 1 | 0.006168583925 | 0.006168502751 | 1.32e−05 |
| 5 | 0.154263309342 | 0.154212568767 | 3.29e−04 |

**Maximum error 3.29e−04 · |Elmer − Python| = 4.21e−11** (machine precision)

### 2D — square well, L=10, 40×40 (1521 free DOFs)

| n | Elmer | analytic | error |
|---|---|---|---|
| 1 | 0.1974935767 | 0.1973920880 | 5.14e−04 |
| 2,3 | 0.4943433728 | 0.4934802201 | 1.75e−03 (degenerate) |
| 4 | 0.7911931689 | 0.7895683521 | 2.06e−03 |
| 5 | 0.9911281920 | 0.9869604401 | 4.22e−03 |

### 3D — cubic well, L=6, 10×10×10 (729 free DOFs)

| n | Elmer | analytic | error |
|---|---|---|---|
| 1 | 0.8292535815 | 0.8224670334 | 8.25e−03 |
| 2,3,4 | 1.6859901662 | 1.6449340668 | 2.50e−02 (triply degenerate) |
| 5 | 2.5427267508 | 2.4674011003 | 3.05e−02 |

### Convergence — order O(h²) confirmed

**2D:**
| n | error | order |
|---:|---|---:|
| 10 | 8.251e−03 | — |
| 20 | 2.058e−03 | 2.00 |
| 40 | 5.141e−04 | 2.00 |
| 80 | 1.285e−04 | 2.00 |

**3D:**
| n | DOFs | error | order |
|---:|---:|---|---:|
| 5 | 64 | 3.331e−02 | — |
| 8 | 343 | 1.292e−02 | 2.16 |
| 10 | 729 | 8.251e−03 | 1.98 |
| 12 | 1331 | 5.725e−03 | 2.04 |
| 15 | 2744 | 3.661e-03 | 1.98 |

---

## 4. Bugs fixed during development

These are errors I introduced myself and had to diagnose — they are worth recording
because they are common pitfalls:

### a) Dirichlet with diagonal = 1 produced spurious eigenvalues
Zeroing the row/column and setting `A(i,i)=B(i,i)=1` creates artificial modes with
**eigenvalue exactly 1.0**. In the 3D case this contaminated the answer (the five
"eigenvalues" came out as 0.829, 1.0, 1.0, 1.0, 1.0).
**Fix:** actual elimination of the prescribed DOFs — extract submatrices
`Af(i,j) = A(IdxFree(i), IdxFree(j))` and solve only the reduced system.

### b) Mesh header format
The second section of `mesh.header` is `element code` + **number of NODES per
element**, not the element count. I mistakenly used `404 160` when the correct
value was `202 2` (a 2-node line).

### c) Boundary element code
2D and 3D boundaries are **lines** with 2 nodes → code **202**. I used `404`
(quadrilateral) by mistake.
3D boundary: quad faces → **404** (correct).

### d) `.sif` with a fixed dimension
A `Coordinate System = Cartesian 2D` applied to a 3D mesh produces
`ERROR:: ElementMetric: Degenerate 3D element` and absurd eigenvalues (−2.8e32).
The dimension in the SIF **must** match that of the mesh.

---

## 5. Delivered tools

| File | Function |
|---|---|
| `fonte/SchrodingerSolver.F90` | the solver (standard Elmer interface) |
| `cases/generate_elmer_mesh.py` | structured 1D/2D/3D mesh generator in Elmer format |
| `verification/run_schrodinger.pbs` | PBS submission on the compute nodes |
| `shares/VERIFICACAO_CRUZADA.txt` | consolidated results |
| `shares/ELMER_caso*.log` | verbatim solver output |

### Why a mesh generator of its own
`ElmerGrid`'s `Minimum Element Divisions` **does not** produce the requested count
(I asked for 39×39 and got 1560 elements instead of 1521). The Python generator gives
exact control and is verified against ElmerGrid's own mesh format.

### Compilation
```bash
module load gnu12 openmpi4 openblas
MODDIR=$HOME/elmer_build/build_gnu/fmodules
LIBDIR=$HOME/elmer_build/build_gnu/fem/src
mpif90 -c SchrodingerSolver.F90 -o SchrodingerSolver.o -fallow-argument-mismatch \
  -DCONTIG= -DHAVE_EXECUTECOMMANDLINE -DUSE_ISO_C_BINDINGS -O2 -fPIC -I$MODDIR
mpif90 -shared SchrodingerSolver.o -o SchrodingerSolver.so -L$LIBDIR -lelmersolver
```

---

## 6. Current limitations

1. **Uniform potential only** — it reads `Potential` from the Solver. Nodal V(x) is
   missing (barrier, double well, sinusoidal potential).
2. **`DSYGV` is serial and dense** — O(n²) memory. The 2D case with 1521 DOFs already
   requires ~40 MB per matrix. For large meshes ARPACK (with the patch) or a parallel
   eigensolver will be needed.
3. **No MPI** — the solver runs in 1 process.
4. The header comment mentions a generic `Eigen System Values`; in practice it
   uses `Apply Dirichlet` + `Eigen System Values`.

---

## 7. Next steps (order defined by the user)

~~2. Test 2D and 3D~~ **COMPLETED** (this round)
→ **1. Variable potential V(x)** (next)
→ 3. MPI parallelization
