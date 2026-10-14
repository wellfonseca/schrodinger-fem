# Raw timing data

The numbers in the performance table of the manuscript, together with a
re-measurement taken in a later session on the same machine. Both are reported
so that the run-to-run variation is visible rather than hidden.

## Method

* 2D harmonic case, `L = 20`, `alpha = 1/4`, `k = 4` eigenpairs requested.
* Sparse route: `sparse_eigsolver.solve_sparse`, ARPACK shift-and-invert with
  `sigma = 0`, `tol = 0`, `ncv` and `maxiter` at library defaults.
* Dense route: `sparse_eigsolver.solve_dense`, LAPACK `DSYGV`.
* Single thread (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`), wall-clock time of
  the eigensolution call alone, excluding assembly.
* Sparse figures are the minimum of three repetitions; dense figures a single
  run, because of their cost.

## Session A (manuscript table)

| n_f | dense (s) | sparse (s) |
|---|---|---|
| 361 | 0.01 | < 0.01 |
| 1 521 | 0.58 | 0.01 |
| 6 241 | 41.06 | 0.05 |
| 14 161 | 491.60 | 0.12 |
| 39 601 | infeasible (37.6 GB) | 0.41 |

## Session B (re-measurement, same machine, later day)

| n_f | dense (s) | sparse min (s) | sparse median (s) | sparse max (s) |
|---|---|---|---|---|
| 361 | 0.03 | 0.009 | 0.012 | 0.013 |
| 1 521 | 1.33 | 0.032 | 0.033 | 0.035 |
| 6 241 | 82.91 | 0.150 | 0.153 | 0.157 |
| 14 161 | 837.90 | 0.376 | 0.378 | 0.382 |
| 39 601 | not attempted | 1.229 | 1.253 | 1.258 |

## Thread scaling (session B, n_f = 39 601, minimum of three)

| OMP_NUM_THREADS | time (s) | ratio |
|---|---|---|
| 1 | 0.993 | 1.00 |
| 2 | 1.029 | 1.04 |
| 4 | 1.273 | 1.28 |
| 8 | 1.901 | 1.91 |

Threading the inner operations gives no speed-up; it degrades the result.

## What this shows

The absolute times differ between sessions by a factor of two to three, with the
machine otherwise idle in both. The speed-up ratios are stable: about 42, 553 and
2 228 in session B against 58, 821 and 4 097 in session A, for the same three
problem sizes. The algorithmic conclusion --- that the sparse route changes the
accessible problem size rather than merely the constant --- is the same in both.
