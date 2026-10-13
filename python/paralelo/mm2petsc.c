/*=============================================================================
 *  mm2petsc.c — converte Matrix Market (coordenado, simetrico/geral)
 *               para o formato binario do PETSc.
 *
 *  Uso:  ./mm2petsc entrada.mtx saida.petsc
 *
 *  Le o cabecalho "%%MatrixMarket matrix coordinate real general|symmetric",
 *  a linha de dimensoes "n n nnz" e as triplas "i j v" (1-based).
 *  Monta a matriz no PETSc e grava em binario paralelo.
 *===========================================================================*/
#include <petsc.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv)
{
  Mat            A;
  PetscViewer    v;
  PetscInt       nrow = 0, ncol = 0, nnz = 0, i, j;
  PetscInt       nloc, rstart;
  PetscScalar    val;
  PetscErrorCode ierr, ierr2;
  FILE          *f;
  char           line[512], object[64], format[64], field[64], symmetry[64];
  PetscMPIInt    rank, size;
  PetscInt       cnt = 0;

  ierr = PetscInitialize(&argc, &argv, NULL, "mm2petsc"); CHKERRQ(ierr);
  if (argc < 3) { PetscPrintf(PETSC_COMM_WORLD, "uso: mm2petsc in.mtx out.petsc\n"); return 1; }
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if (size > 1) {
    PetscPrintf(PETSC_COMM_WORLD, "aviso: usar 1 processo para converter\n");
  }

  f = fopen(argv[1], "r");
  if (!f) { PetscPrintf(PETSC_COMM_WORLD, "erro: nao abriu %s\n", argv[1]); return 1; }

  /* cabecalho */
  if (!fgets(line, sizeof(line), f)) { fclose(f); return 1; }
  if (sscanf(line, "%%%%MatrixMarket %63s %63s %63s %63s",
             object, format, field, symmetry) < 4) {
    PetscPrintf(PETSC_COMM_WORLD, "erro: cabecalho invalido\n"); fclose(f); return 1;
  }
  /* comentarios */
  while (fgets(line, sizeof(line), f) && line[0] == '%') { }
  sscanf(line, "%d %d %d", &nrow, &ncol, &nnz);
  PetscPrintf(PETSC_COMM_WORLD, "  lido: %d x %d, nnz = %d, simetria = %s\n",
              (int)nrow, (int)ncol, (int)nnz, symmetry);

  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol); CHKERRQ(ierr);
  ierr = MatSetType(A, MATAIJ); CHKERRQ(ierr);
  ierr = MatSetUp(A); CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A, &rstart, NULL); CHKERRQ(ierr);
  ierr = MatGetLocalSize(A, &nloc, NULL); CHKERRQ(ierr);

  while (fgets(line, sizeof(line), f)) {
    int a, b; double c;
    if (line[0] == '%' || line[0] == '\n') continue;
    if (sscanf(line, "%d %d %lf", &a, &b, &c) != 3) continue;
    i = a - 1; j = b - 1; val = c;
    if (i >= rstart && i < rstart + nloc) {
      ierr = MatSetValue(A, i, j, val, INSERT_VALUES); CHKERRQ(ierr);
      if (strncasecmp(symmetry, "symmetric", 9) == 0 && i != j) {
        ierr = MatSetValue(A, j, i, val, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
    cnt++;
  }
  fclose(f);
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "  entradas lidas: %d\n", (int)cnt);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, argv[2], FILE_MODE_WRITE, &v); CHKERRQ(ierr);
  ierr = MatView(A, v); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "  gravado: %s\n", argv[2]);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
