# Potential syntax in Elmer — and the multidimensional MATC obstacle

**Reference document.** It records the correct syntax for defining `V(x,y,z)`
in the SIF, including the one that **takes time to discover**.

---

## 1. The essentials

### 1D — works directly
```
Body Force 1
  Potential = Variable Coordinate 1
    Real MATC "1.0*tx*tx"
End
```
With **one** variable, `tx` is its value.

### 2D / 3D — USE A COMMA as separator
```
Body Force 1
  Potential = Variable Coordinate 1, Coordinate 2, Coordinate 3
    Real MATC "0.25*(tx*tx + tx(1)*tx(1) + tx(2)*tx(2))"
End
```
With **several** variables, `tx` is the **first** one and `tx(1)`, `tx(2)`, … are the
following ones. (Manual: *"We could have more than one variable. In that case, we use
the vector notation (example: tx(1) or tx(2))"*.)

**Verified:** 3D harmonic, box L=8, mesh 14³
→ `V_min = 0`, `V_max = 12.0` = `0.25·(16+16+16)` ✅
→ `E₀ = 1.558551`, **triple** degeneracy at `E₁ = 2.591398` ✅

---

## 2. The syntaxes that do NOT work (all tested)

| Syntax | Result |
|---|---|
| `Variable Coordinate 1 2` | `Can't find independent variable:[coordinate 1 2]` |
| `Variable Coordinate 1` + `tx(1)` | `MATC ERROR: Index out of bounds` |
| `Variable Coordinate 1 Coordinate 2` (space) | `Can't find independent variable` |
| two keywords `Xc`/`Yc` with distinct names | `MATC ERROR: Undeclared identifier` |
| multi-line `Real MATC` (3 lines) | accepted but **with no effect** |

### Why the space fails — the root cause

In the Elmer source (`Lists.F90`, `ListParseStrToVars`), the separator is a
**COMMA**, not a space:

```fortran
l1 = INDEX( str(l0:slen), ',' )     ! <-- separator is a COMMA
```

Without the comma, `Coordinate 1 Coordinate 2` is treated as a **single variable
name** — which does not exist.

> **Lesson recorded:** the manual shows `Velocity 1 = Variable Coordinate 2`
> (full name) but **does not make the separator clear**. The answer was in the
> source code, not in the documentation.

---

## 3. Alternative that also works (implemented in the solver)

The solver also accepts **component-wise** construction, without multidimensional MATC:

```
Body Force 1
  Potential X = Variable Coordinate 1
    Real MATC "0.25*tx*tx"
  Potential Y = Variable Coordinate 2
    Real MATC "0.25*tx*tx"
  Potential Z = Variable Coordinate 3
    Real MATC "0.25*tx*tx"
End
```
The solver adds `V = Vx + Vy + Vz` internally. Useful because each term is 1D and
can come from different sources.

**Precedence order in the solver:**
1. `Potential` (nodal field — accepts both syntaxes above)
2. `Potential X` + `Potential Y` + `Potential Z`
3. `V = 0`

---

## 4. Where this is in the source

| File | What it contains |
|---|---|
| `Lists.F90` → `ListParseStrToVars` | parser for `Variable`, separator `,` |
| `Lists.F90` → `ListGetReal` | evaluates the expression per node |
| `DefUtils.F90` → `GetRealValues` | element-wise variant (do not use outside the loop) |
| `EigenSolve.F90` → `ArpackEigenSolve` | **eigensolver with a bug** (`RESID` not initialized) |

Reference copies of the Elmer 9.0 source:
`<ElmerRun>/...`
(reference **read-only** — local copy in `referencias/elmer/`)

---

## 5. Credit

The reference source is in `<Elmer installation>` (Elmer 9.0 with MUMPS,
installed locally) — **read-only**. The `ElmerRun/trafo2026/` folder
already used the `.elmer_ref/` pattern for API consultation, which I adopted here.
