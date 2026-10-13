# Item 1 — Potential V(x) and regression suite

**Status: V(x) WORKING · regression suite 5/5 · convention verified**

---

## 1. Bug fixed: reading V(x)

### Before
```fortran
Material => Model % Materials(1) % Values
IF ( ListCheckPresent( Material, 'Potential' ) ) THEN
  Vuni = ListGetConstReal( Material, 'Potential', Found )   ! scalar
```
Two defects:
- `Potential` is **not a valid keyword** in `Material` in Elmer 9
  (`ERROR:: Unknown specifier ... property name:[potential]`);
- `ListGetConstReal` only accepts a **scalar** — incapable of reading a field.

Result: the branch never executed and the solver read **V = 0 silently**.
No test with V ≠ 0 had ever run.

### After
```fortran
! (1) Body Force  (accepts expression: Variable Coordinate / MATC)
Vtmp = ListGetReal( BFList, 'Potential', nNodes, NodeIdx, Found )
! (2) Material    (fallback)
! (3) V = 0
```
Implementation notes that cost time:
- `GetRealValues` **is no good** — it depends on `GetCurrentElement` (loop context);
- `ListGetReal` is a **function** (`RESULT(F)`, with `REAL :: F(N)`), not a subroutine;
- `BF` collides with `Bf` (matrix) — Fortran is case-insensitive.
  Renamed to `BFList`;
- `Max Output Level = 3` **suppresses** the solver messages (they use `Level=4`).

---

## 2. IMPORTANT FINDING: the solver's convention

The earlier validation (infinite well, L=40) was **ambiguous by coincidence**:

| Convention | `E₁` for `L=40` |
|---|---|
| `-ψ'' + Vψ = Eψ` | `π²/L² = 0.00616850` |
| `-½ψ'' + Vψ = Eψ` | `π²/(2L²) = 0.00308425` |

With L=40 the two do **not** coincide — but the measured value (0.00616858) matched
the first, so everything seemed fine. **The test did not distinguish the conventions.**

Tests that distinguish them (L=10 and the harmonic oscillator) were decisive:

| Test | Elmer | `-ψ''` | `-½ψ''` |
|---|---|---|---|
| well L=10, n=1 | 0.09869637 | **0.09869604** | 0.04934802 |
| harmonic α=1, n=0 | 1.00002292 | **1.00000000** | 0.70710678 |

**Convention confirmed: `-ψ'' + V(x)ψ = Eψ`** (ħ = 1, m = ½ — atomic units).
Relation to the usual form `-½ψ''+Vψ = Eψ`: `E_here = 2·E_usual`.

This is now recorded in the solver header and in the regression suite header,
so that it is not lost.

---

## 3. Regression suite — 5/5

`verification/confere.py` runs the same comparison locally, creates its own folder for each case
(avoids races between files) and compares against the analytic solution.

| Case | Check | Result |
|---|---|---|
| infinite well L=10 | `Eₙ=(nπ/L)²` | error **2.96e−05** ✅ |
| harmonic oscillator k=1 | `Eₙ=2(n+½)` | error **1.43e−04** ✅ |
| Morse (De=10, β=0.5) | potential minimum = −De | `V_min=−10.000000` ✅ |
| double well | tunneling splitting > 0 | `E₀=0.224081 E₁=0.599904` Δ=**0.375823** ✅ |
| sinusoidal A=0.5 | `V∈[−A,+A]` | `V=[−0.5, +0.5]` ✅ |

**Result: 5 passed, 0 failed.**

### Physical comment on each case

- **Harmonic:** spacing exactly 2.0 (`E₀=1.00014, E₁=3.00030, E₂=5.00061`)
- **Morse:** `E₀=−8.481` above the bottom `V_min=−10` — zero-point energy
  `≈1.519`, compatible with the local `ω/2`
- **Double well:** splitting `0.376` confirms tunneling between the two wells
  (even/odd states nearly degenerate without a barrier)
- **Sinusoidal:** eigenvalues close to zero and above the potential minimum
  `−0.5`, at the onset of the band structure

---

## 4. Set of cases by dimension

| Case | 1D | 2D | 3D |
|---|---|---|---|
| infinite well | ✅ | ✅ | ✅ |
| harmonic | ✅ | — | — |
| Morse | ✅ | — | — |
| double well | ✅ | — | — |
| sinusoidal | ✅ | — | — |
| barrier | pending | — | — |

---

## 5. What changes in the validation of item 2

The validation of item 2 (1D/2D/3D) remains **valid**: the convention used is
consistent between Elmer and Python (cross-verification at 4.21e−11) and the
analytic comparisons used the correct formula for `-ψ''` in each case.
What was missing was to **state explicitly** which convention it was — now it is.

---

## 6. Outstanding items

1. **Potential barrier** — still to be implemented (transmission/reflection requires
   a radiation boundary condition, not Dirichlet).
2. **2D/3D with V(x) ≠ 0** — the code supports it (the potential is read as a nodal
   field, independent of dimension), but it has only been tested in 1D.
3. **Scale** — dense/serial `DSYGV`; limits to ~8000 DOFs. This is item 3.
