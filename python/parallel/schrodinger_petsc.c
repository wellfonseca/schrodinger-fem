/*=============================================================================
 *  schrodinger_petsc.c
 *
 *  PARALLEL eigensolver for the generalised problem  A x = E M x
 *  (stationary Schrodinger equation by finite elements).
 *
 *  Reads the matrices A and M assembled by the FEM (PETSc binary format),
 *  solves for the k smallest eigenvalues with SLEPc in shift-invert mode
 *  (Krylov-Schur + PETSc parallel LU factorization) and prints the results.
 *
 *  Compile:
 *    module load gnu12 openmpi4 petsc slepc
 *    mpicc schrodinger_petsc.c -o schrodinger_petsc \
 *          -I$PETSC_DIR/include -I$SLEPC_DIR/include \
 *          -L$PETSC_DIR/lib -L$SLEPC_DIR/lib -lslepc -lpetsc -lm
 *
 *  Run:
 *    mpirun -np N ./schrodinger_petsc -A A.petsc -M M.petsc -k 4
 *
 *  Output: eigenvalues + timing metrics per stage.
 *===========================================================================*/
#include <petsc.h>
#include <slepc.h>

static char help[] = "Parallel FEM eigensolver: A x = E M x\n";

int main(int argc, char **argv)
{
  Mat            A, M;
  EPS            eps;
  ST             st;
  KSP            ksp;
  PC             pc;
  Vec            xr, xi;
  PetscInt       k = 4, nconv, i;
  PetscScalar    kr, ki;
  PetscReal      err;
  PetscErrorCode ierr;
  char           fileA[PETSC_MAX_PATH_LEN] = "A.petsc";
  char           fileM[PETSC_MAX_PATH_LEN] = "M.petsc";
  PetscLogDouble t0, t1, t2, t3;
  PetscMPIInt    size, rank;

  ierr = SlepcInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  ierr = PetscOptionsGetString(NULL, NULL, "-A", fileA, sizeof(fileA), NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, "-M", fileM, sizeof(fileM), NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-k", &k, NULL); CHKERRQ(ierr);

  ierr = PetscTime(&t0); CHKERRQ(ierr);

  /* --- read the matrices --- */
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetType(A, MATAIJ); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A, "A"); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  {
    PetscViewer v;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileA, FILE_MODE_READ, &v); CHKERRQ(ierr);
    ierr = MatLoad(A, v); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
  }

  ierr = MatCreate(PETSC_COMM_WORLD, &M); CHKERRQ(ierr);
  ierr = MatSetType(M, MATAIJ); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)M, "M"); CHKERRQ(ierr);
  ierr = MatSetFromOptions(M); CHKERRQ(ierr);
  {
    PetscViewer v;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileM, FILE_MODE_READ, &v); CHKERRQ(ierr);
    ierr = MatLoad(M, v); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
  }

  if (rank == 0) {
    PetscPrintf(PETSC_COMM_WORLD,
      "==================================================\n");
    PetscPrintf(PETSC_COMM_WORLD,
      " Parallel Schrodinger FEM (SLEPc shift-invert)\n");
    PetscPrintf(PETSC_COMM_WORLD,
      " MPI processes : %d\n", (int)size);
  }
  {
    MatInfo     info;
    PetscInt    nrow, ncol;
    ierr = MatGetInfo(A, MAT_GLOBAL_SUM, &info); CHKERRQ(ierr);
    ierr = MatGetSize(A, &nrow, &ncol); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
        " matrix A      : %d x %d,  nnz = %.0f\n",
        (int)nrow, (int)ncol, (double)info.nz_used); CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " eigenvalues k : %d\n", (int)k); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "==================================================\n"); CHKERRQ(ierr);

  ierr = PetscTime(&t1); CHKERRQ(ierr);

  /* --- eigensolver --- */
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, A, M); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_GHEP); CHKERRQ(ierr);
  ierr = EPSSetType(eps, EPSKRYLOVSCHUR); CHKERRQ(ierr);
  /* target-aware: required for shift-and-invert in SLEPc */
  ierr = EPSSetTarget(eps, 0.0); CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, k, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps, 1.0e-10, 500); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  /* shift-invert: parallel LU factorization and reuse */
  ierr = EPSGetST(eps, &st); CHKERRQ(ierr);
  ierr = STSetType(st, STSINVERT); CHKERRQ(ierr);
  ierr = STGetKSP(st, &ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  /* PETSc LU factorizer.
     Not every build ships MUMPS/SuperLU/UMFPACK. On this cluster only the
     internal ("petsc") LU is available -- verified with MatGetFactor.
     It can be switched from the command line:
        -fact_solver mumps|superlu_dist|superlu|umfpack|petsc
     We do NOT probe automatically: MatGetFactor emits the PETSc error
     before returning the code, which aborts execution. */
  {
    char fact[64] = "petsc";
    PetscBool flg;
    ierr = PetscOptionsGetString(NULL, NULL, "-fact_solver", fact, sizeof(fact), &flg);
    CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverType(pc, fact); CHKERRQ(ierr);
    if (rank == 0)
      ierr = PetscPrintf(PETSC_COMM_WORLD, " LU factorizer : %s\n", fact); CHKERRQ(ierr);
  }
  ierr = STSetShift(st, 0.0); CHKERRQ(ierr);

  ierr = EPSSolve(eps); CHKERRQ(ierr);
  ierr = PetscTime(&t2); CHKERRQ(ierr);

  ierr = EPSGetConverged(eps, &nconv); CHKERRQ(ierr);
  ierr = MatCreateVecs(A, &xr, &xi); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "\n--- EIGENVALUES (energies) ---\n"); CHKERRQ(ierr);
  for (i = 0; i < nconv; i++) {
    ierr = EPSGetEigenpair(eps, i, &kr, &ki, xr, xi); CHKERRQ(ierr);
    ierr = EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &err); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
        "  E(%d) = %+.12e   (rel. error %.2e)\n",
        (int)(i+1), (double)PetscRealPart(kr), (double)err); CHKERRQ(ierr);
  }
  ierr = PetscTime(&t3); CHKERRQ(ierr);

  /* --- metrics --- */
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "\n--- PERFORMANCE ---\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " read matrices    : %8.3f s\n", (double)(t1-t0)); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " SLEPc solve      : %8.3f s\n", (double)(t2-t1)); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " extraction/eval  : %8.3f s\n", (double)(t3-t2)); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " TOTAL            : %8.3f s   with %d processes\n",
      (double)(t3-t0), (int)size); CHKERRQ(ierr);

  VecDestroy(&xr); VecDestroy(&xi);
  EPSDestroy(&eps);
  MatDestroy(&A); MatDestroy(&M);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
