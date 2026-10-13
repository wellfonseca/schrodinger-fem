/*=============================================================================
 *  schrodinger_petsc.c
 *
 *  Autossolver PARALELO para o problema generalizado  A x = E M x
 *  (Schrodinger estacionario por elementos finitos).
 *
 *  Le as matrizes A e M montadas pelo FEM (formato PETSc binario),
 *  resolve os k menores autovalores com SLEPc em modo shift-invert
 *  (Krylov-Schur + fatoracao LU paralela do PETSc) e imprime os resultados.
 *
 *  Compilar:
 *    module load gnu12 openmpi4 petsc slepc
 *    mpicc schrodinger_petsc.c -o schrodinger_petsc \
 *          -I$PETSC_DIR/include -I$SLEPC_DIR/include \
 *          -L$PETSC_DIR/lib -L$SLEPC_DIR/lib -lslepc -lpetsc -lm
 *
 *  Executar:
 *    mpirun -np N ./schrodinger_petsc -A A.petsc -M M.petsc -k 4
 *
 *  Saida: autovalores + metricas de tempo por etapa.
 *===========================================================================*/
#include <petsc.h>
#include <slepc.h>

static char help[] = "Autossolver paralelo FEM: A x = E M x\n";

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

  /* --- le as matrizes --- */
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
      " Schrodinger FEM paralelo (SLEPc shift-invert)\n");
    PetscPrintf(PETSC_COMM_WORLD,
      " processos MPI : %d\n", (int)size);
  }
  {
    MatInfo     info;
    PetscInt    nrow, ncol;
    ierr = MatGetInfo(A, MAT_GLOBAL_SUM, &info); CHKERRQ(ierr);
    ierr = MatGetSize(A, &nrow, &ncol); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
        " matriz A      : %d x %d,  nnz = %.0f\n",
        (int)nrow, (int)ncol, (double)info.nz_used); CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " autovalores k : %d\n", (int)k); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "==================================================\n"); CHKERRQ(ierr);

  ierr = PetscTime(&t1); CHKERRQ(ierr);

  /* --- autossolver --- */
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, A, M); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_GHEP); CHKERRQ(ierr);
  ierr = EPSSetType(eps, EPSKRYLOVSCHUR); CHKERRQ(ierr);
  /* ciente do alvo: obrigatorio para shift-and-invert no SLEPc */
  ierr = EPSSetTarget(eps, 0.0); CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, k, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps, 1.0e-10, 500); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  /* shift-invert: fatoracao LU paralela e reuso */
  ierr = EPSGetST(eps, &st); CHKERRQ(ierr);
  ierr = STSetType(st, STSINVERT); CHKERRQ(ierr);
  ierr = STGetKSP(st, &ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  /* Fatorador LU do PETSc.
     Nem toda build traz MUMPS/SuperLU/UMFPACK. Neste cluster apenas o LU
     interno ("petsc") esta disponivel -- verificado com MatGetFactor.
     Pode ser trocado por linha de comando:
        -fact_solver mumps|superlu_dist|superlu|umfpack|petsc
     NAO fazemos probing automatico: MatGetFactor emite o erro do PETSc
     antes de devolver o codigo, o que aborta a execucao. */
  {
    char fact[64] = "petsc";
    PetscBool flg;
    ierr = PetscOptionsGetString(NULL, NULL, "-fact_solver", fact, sizeof(fact), &flg);
    CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverType(pc, fact); CHKERRQ(ierr);
    if (rank == 0)
      ierr = PetscPrintf(PETSC_COMM_WORLD, " fatorador LU  : %s\n", fact); CHKERRQ(ierr);
  }
  ierr = STSetShift(st, 0.0); CHKERRQ(ierr);

  ierr = EPSSolve(eps); CHKERRQ(ierr);
  ierr = PetscTime(&t2); CHKERRQ(ierr);

  ierr = EPSGetConverged(eps, &nconv); CHKERRQ(ierr);
  ierr = MatCreateVecs(A, &xr, &xi); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "\n--- AUTOVALORES (energias) ---\n"); CHKERRQ(ierr);
  for (i = 0; i < nconv; i++) {
    ierr = EPSGetEigenpair(eps, i, &kr, &ki, xr, xi); CHKERRQ(ierr);
    ierr = EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &err); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
        "  E(%d) = %+.12e   (erro rel. %.2e)\n",
        (int)(i+1), (double)PetscRealPart(kr), (double)err); CHKERRQ(ierr);
  }
  ierr = PetscTime(&t3); CHKERRQ(ierr);

  /* --- metricas --- */
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "\n--- DESEMPENHO ---\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " leitura matrizes : %8.3f s\n", (double)(t1-t0)); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " solucao SLEPc    : %8.3f s\n", (double)(t2-t1)); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " extracao/eval    : %8.3f s\n", (double)(t3-t2)); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      " TOTAL            : %8.3f s   com %d processos\n",
      (double)(t3-t0), (int)size); CHKERRQ(ierr);

  VecDestroy(&xr); VecDestroy(&xi);
  EPSDestroy(&eps);
  MatDestroy(&A); MatDestroy(&M);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
