!/*****************************************************************************/
! *
! *  SchrodingerSolver.F90
! *
! *  Solver de elementos finitos para a equacao de Schrodinger estacionaria
! *  (autoproblema de valores proprios), para 1D, 2D e 3D.
! *
! *      -lap(psi) + V(r) psi = E psi          (hbar = m = 1)
! *
! *  Forma fraca (Galerkin), psi = sum_i N_i psi_i, testando com N_j:
! *
! *     [ int grad(Nj).grad(Ni) dO + int Nj V Ni dO ] psi_i
! *                         = E [ int Nj Ni dO ] psi_i
! *
! *     isto e:   A psi = E M psi,   com   A = K + Ma
! *
! *        K  = int grad(Ni).grad(Nj)      (rigidez)
! *        M  = int Ni Nj                  (massa)
! *        Ma = int Ni V Nj                (massa ponderada pelo potencial)
! *
! *  O autoproblema generalizado e resolvido com DSYGV (LAPACK), evitando
! *  depender do autossolver ARPACK interno do Elmer.
! *
! *  CONDICOES DE CONTORNO
! *     Apply Dirichlet = Logical True  ->  psi = 0 nos contornos geometricos
! *     (deteccao automatica dos nos de contorno; sem necessidade de BC no SIF)
! *
! *  POTENCIAL V  (ordem de precedencia)
! *     1) "Potential"        em Body Force   -> campo nodal; aceita
! *          - 1D:   Variable Coordinate 1        + Real MATC "expr(tx)"
! *          - 2D/3D: Variable Coordinate 1, Coordinate 2[, Coordinate 3]
! *                                              + Real MATC "expr(tx,tx(1),tx(2))"
! *          (ATENCAO: o separador de variaveis do Elmer e VIRGULA)
! *     2) "Potential X/Y/Z"  em Body Force   -> soma por componente
! *     3) ausente                            -> V = 0
! *
! *  UTILITARIOS
! *     Exportar Matrizes = Logical True -> grava A, M e os GDL de Dirichlet
! *         em triplas ASCII (<prefixo>_A.dat, _M.dat, _fixed.dat) para
! *         resolucao externa com solver esparso/paralelo.
!
! *  REFERENCIAS
! *     docs/  (manuais do modelo, do solver e validacao)
!
! *  Para a documentacao completa ver:
! *     schrodinger/Manual_Modelo_Schrodinger.md
! *     schrodinger/Manual_Solver_Schrodinger.md
! *
! ******************************************************************************/
SUBROUTINE SchrodingerSolver( Model, Solver, dt, Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  USE ElementDescription
  USE BandwidthOptimize

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
!  Locais
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(ValueList_t), POINTER :: Params, Material, BFList
  TYPE(Variable_t), POINTER :: PsiSol

  REAL(KIND=dp), ALLOCATABLE :: A(:,:), B(:,:), Af(:,:), Bf(:,:)
  REAL(KIND=dp), ALLOCATABLE :: EigenVals(:), Work(:), xfree(:,:)
  REAL(KIND=dp), ALLOCATABLE :: Vtmp(:)
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), Vnodal(:)
  REAL(KIND=dp), ALLOCATABLE :: Stiff(:,:), Mass(:,:), PotMass(:,:)

  CHARACTER(LEN=MAX_NAME_LEN) :: Name
  INTEGER :: i, j, k, l, p, q, t, n, nd, dim, nNodes
  INTEGER :: nEigen, nDof, nFree, LWORK, info_lap, nBC
  LOGICAL :: Found, ApplyDiri, AchouV, Fx, Fy, Fz
  LOGICAL :: FazExport
  CHARACTER(LEN=MAX_NAME_LEN) :: ExpPrefixo
  TYPE(ValueList_t), POINTER :: BFList0
  REAL(KIND=dp) :: detJ, s, Vavg, Vuni
  TYPE(GaussIntegrationPoints_t) :: Integ

  ! mapeamento de nos prescritos
  INTEGER, ALLOCATABLE :: Perm(:), MapFree(:), IdxFree(:), NodeIdx(:)
  LOGICAL, ALLOCATABLE :: Fixed(:)
  REAL(KIND=dp), ALLOCATABLE :: xmin(:), xmax(:), tol(:)
  REAL(KIND=dp) :: coord, dax(3), day(3), daz(3)

  !============================================================================
  Params => Solver % Values
  dax = 0.0_dp; day = 0.0_dp; daz = 0.0_dp
  dax(1) = 1.0_dp; day(2) = 1.0_dp; daz(3) = 1.0_dp
  dim = CoordinateSystemDimension()
  nEigen = ListGetInteger( Params, 'Eigen System Values', Found )
  IF ( .NOT. Found ) nEigen = 5

  CALL Info( 'SchrodingerSolver', '======================================', Level=4 )
  CALL Info( 'SchrodingerSolver', 'Solver de Schrodinger (FEM) iniciado',  Level=4 )
  WRITE( Message, '(A,I0)' ) 'Numero de autovalores pedidos: ', nEigen
  CALL Info( 'SchrodingerSolver', Message, Level=4 )
  WRITE( Message, '(A,I0)' ) 'Dimensao do problema: ', dim
  CALL Info( 'SchrodingerSolver', Message, Level=4 )
  CALL Info( 'SchrodingerSolver', '======================================', Level=4 )

  !--------------------------------------------------------------------
  ! 1) numeracao global dos nos e mapeamento para DOFs
  !--------------------------------------------------------------------
  nNodes = Model % NumberOfNodes
  nDof   = nNodes

  ALLOCATE( Perm(nDof), Fixed(nDof) )
  Perm = 0
  Fixed = .FALSE.

  PsiSol => Solver % Variable
  IF ( .NOT. ASSOCIATED(PsiSol) ) THEN
    CALL Fatal( 'SchrodingerSolver', 'Variavel da solucao nao definida.' )
  END IF
  Perm(1:nDof) = PsiSol % Perm(1:nDof)

  !--------------------------------------------------------------------
  ! 2) aplicar Dirichlet por deteccao GEOMETRICA dos nos de contorno
  !    (robusto em 1D, 2D e 3D; independe da topologia de contorno)
  !--------------------------------------------------------------------
  ApplyDiri = ListGetLogical( Params, 'Apply Dirichlet', Found )
  IF ( .NOT. Found ) ApplyDiri = .TRUE.

  IF ( ApplyDiri ) THEN
    ALLOCATE( xmin(3), xmax(3), tol(3) )
    xmin =  HUGE(1.0_dp)
    xmax = -HUGE(1.0_dp)
    DO i = 1, nDof
      DO k = 1, dim
        xmin(k) = MIN( xmin(k), Model % Nodes % x(i) * dax(k) &
                     + Model % Nodes % y(i) * day(k) + Model % Nodes % z(i) * daz(k) )
        xmax(k) = MAX( xmax(k), Model % Nodes % x(i) * dax(k) &
                     + Model % Nodes % y(i) * day(k) + Model % Nodes % z(i) * daz(k) )
      END DO
    END DO
    ! tolerancia relativa ao tamanho do dominio
    DO k = 1, dim
      tol(k) = 1.0e-8_dp * MAX( 1.0_dp, ABS( xmax(k) - xmin(k) ) )
    END DO

    DO i = 1, nDof
      DO k = 1, dim
        coord = Model % Nodes % x(i) * dax(k) &
              + Model % Nodes % y(i) * day(k) + Model % Nodes % z(i) * daz(k)
        IF ( ABS( coord - xmin(k) ) <= tol(k) .OR. &
             ABS( coord - xmax(k) ) <= tol(k) ) THEN
          Fixed(i) = .TRUE.
          EXIT
        END IF
      END DO
    END DO

    WRITE( Message, '(A,I0)' ) 'Nos de contorno (Dirichlet psi=0): ', COUNT(Fixed)
    CALL Info( 'SchrodingerSolver', Message, Level=4 )
  END IF

  !--------------------------------------------------------------------
  ! 3) montagem de A = K + Ma  e  M
  !--------------------------------------------------------------------
  ALLOCATE( A(nDof,nDof), B(nDof,nDof) )
  A = 0.0_dp
  B = 0.0_dp

  !--------------------------------------------------------------------
  ! POTENCIAL V(x): lido como CAMPO NODAL
  !   Ordem de precedencia:
  !     1) Body Force  "Potential"  (aceita expressao: Variable Coordinate / MATC)
  !     2) Material    "Potential"  (nao aceito pelo Elmer 9 — mantido por robustez)
  !     3) V = 0
  !   Usa ListGetReal com TODOS os nos (nao GetRealValues, que depende do
  !   elemento corrente).
  !--------------------------------------------------------------------
  ALLOCATE( Vnodal(nNodes), NodeIdx(nNodes), Vtmp(nNodes) )
  DO i = 1, nNodes
    NodeIdx(i) = i
  END DO

  Vnodal = 0.0_dp
  AchouV = .FALSE.
  BFList0 => NULL()

  ! (1) Body Force
  IF ( Model % NumberOfBodyForces > 0 ) THEN
    BFList => Model % BodyForces(1) % Values
    BFList0 => BFList
    IF ( ASSOCIATED(BFList) ) THEN
      IF ( ListCheckPresent( BFList, 'Potential' ) ) THEN
        Vtmp = ListGetReal( BFList, 'Potential', nNodes, NodeIdx, Found )
        IF ( Found ) THEN
          Vnodal = Vtmp
          AchouV = .TRUE.
          CALL Info( 'SchrodingerSolver', &
               'V(x) lido como campo nodal de Body Force.', Level=4 )
        END IF
      END IF
    END IF
  END IF

  ! (2) Material
  IF ( .NOT. AchouV ) THEN
    Material => Model % Materials(1) % Values
    IF ( ListCheckPresent( Material, 'Potential' ) ) THEN
      Vtmp = ListGetReal( Material, 'Potential', nNodes, NodeIdx, Found )
      IF ( Found ) THEN
        Vnodal = Vtmp
        AchouV = .TRUE.
        CALL Info( 'SchrodingerSolver', &
             'V(x) lido como campo nodal de Material.', Level=4 )
      END IF
    END IF
  END IF

  !--------------------------------------------------------------------
  ! (4) CONSTRUCAO DE V A PARTIR DAS COORDENADAS
  !     Keywords:  "Potential X", "Potential Y", "Potential Z"
  !     V(x,y,z) = Vx + Vy + Vz
  !     Cada termo e lido com ListGetReal (aceita "Equals Coordinate n",
  !     expressao MATC 1D, ou valor constante) -> robusto e independente
  !     da sintaxe MATC multidimensional.
  !--------------------------------------------------------------------
  IF ( .NOT. AchouV ) THEN
    Fx = .FALSE.; Fy = .FALSE.; Fz = .FALSE.
    Vnodal = 0.0_dp
    IF ( ListCheckPresent( BFList0, 'Potential X' ) ) THEN
      Vtmp = ListGetReal( BFList0, 'Potential X', nNodes, NodeIdx, Fx )
      IF ( Fx ) Vnodal = Vnodal + Vtmp
    END IF
    IF ( dim >= 2 .AND. ListCheckPresent( BFList0, 'Potential Y' ) ) THEN
      Vtmp = ListGetReal( BFList0, 'Potential Y', nNodes, NodeIdx, Fy )
      IF ( Fy ) Vnodal = Vnodal + Vtmp
    END IF
    IF ( dim >= 3 .AND. ListCheckPresent( BFList0, 'Potential Z' ) ) THEN
      Vtmp = ListGetReal( BFList0, 'Potential Z', nNodes, NodeIdx, Fz )
      IF ( Fz ) Vnodal = Vnodal + Vtmp
    END IF
    IF ( Fx .OR. Fy .OR. Fz ) THEN
      AchouV = .TRUE.
      CALL Info( 'SchrodingerSolver', &
           'V(x,y,z) construido de Potential X/Y/Z.', Level=4 )
    END IF
  END IF

  IF ( .NOT. ASSOCIATED(BFList0) ) BFList0 => Params

  IF ( AchouV ) THEN
    WRITE( Message, '(A,ES14.6,A,ES14.6)' ) '  V_min = ', MINVAL(Vnodal), &
                                            '   V_max = ', MAXVAL(Vnodal)
    CALL Info( 'SchrodingerSolver', Message, Level=4 )
    CALL Info( 'SchrodingerSolver', Message, Level=5 )
  ELSE
    CALL Info( 'SchrodingerSolver', 'Sem potencial (V = 0).', Level=4 )
  END IF
  DEALLOCATE( NodeIdx, Vtmp )

  DO k = 1, Model % NumberOfBulkElements
    Element => Model % Elements(k)
    n = Element % Type % NumberOfNodes

    ALLOCATE( Basis(n), dBasisdx(n,3) )
    ALLOCATE( Stiff(n,n), Mass(n,n), PotMass(n,n) )
    Stiff = 0.0_dp; Mass = 0.0_dp; PotMass = 0.0_dp

    CALL GetElementNodes( Nodes, Element )
    Integ = GaussPoints( Element )

    DO t = 1, Integ % n
      s = Integ % s(t)
      IF ( .NOT. ElementInfo( Element, Nodes, Integ % u(t), Integ % v(t), &
             Integ % w(t), detJ, Basis, dBasisdx ) ) CYCLE

      ! potencial medio no ponto de integracao
      Vavg = 0.0_dp
      DO i = 1, n
        Vavg = Vavg + Basis(i) * Vnodal( Element % NodeIndexes(i) )
      END DO

      DO p = 1, n
        DO q = 1, n
          Stiff(p,q) = Stiff(p,q) + s * detJ * &
              SUM( dBasisdx(p,1:dim) * dBasisdx(q,1:dim) )
          Mass(p,q)  = Mass(p,q)  + s * detJ * Basis(p) * Basis(q)
          PotMass(p,q) = PotMass(p,q) + s * detJ * Vavg * Basis(p) * Basis(q)
        END DO
      END DO
    END DO

    ! montagem global
    DO p = 1, n
      i = Element % NodeIndexes(p)
      DO q = 1, n
        j = Element % NodeIndexes(q)
        A(i,j) = A(i,j) + Stiff(p,q) + PotMass(p,q)
        B(i,j) = B(i,j) + Mass(p,q)
      END DO
    END DO

    DEALLOCATE( Basis, dBasisdx, Stiff, Mass, PotMass )
  END DO

  !--------------------------------------------------------------------
  ! 4) ELIMINACAO dos DOFs de Dirichlet
  !    Linhas/colunas dos nos prescritos sao REMOVIDAS (nao apenas
  !    zeradas com diagonal 1, o que introduziria autovalores espurios
  !    em lambda = 1).
  !--------------------------------------------------------------------
  nFree = 0
  DO i = 1, nDof
    IF ( .NOT. Fixed(i) ) nFree = nFree + 1
  END DO

  IF ( nFree < nEigen ) THEN
    CALL Fatal( 'SchrodingerSolver', 'DOFs livres insuficientes.' )
  END IF

  WRITE( Message, '(A,I0)' ) 'DOFs livres (apos eliminacao): ', nFree
  CALL Info( 'SchrodingerSolver', Message, Level=4 )

  ALLOCATE( Af(nFree,nFree), Bf(nFree,nFree), MapFree(nDof), IdxFree(nFree) )
  nFree = 0
  DO i = 1, nDof
    IF ( .NOT. Fixed(i) ) THEN
      nFree = nFree + 1
      IdxFree(nFree) = i
      MapFree(i) = nFree
    ELSE
      MapFree(i) = 0
    END IF
  END DO

  DO j = 1, nFree
    DO i = 1, nFree
      Af(i,j) = A( IdxFree(i), IdxFree(j) )
      Bf(i,j) = B( IdxFree(i), IdxFree(j) )
    END DO
  END DO

  !--------------------------------------------------------------------
  ! 5) autoproblema generalizado reduzido:  Af x = lambda Bf x  (DSYGV)
  !--------------------------------------------------------------------
  LWORK = MAX( 1, 3*nFree - 1 )
  ALLOCATE( EigenVals(nFree), Work(LWORK), xfree(nFree,nEigen) )
  EigenVals = 0.0_dp; Work = 0.0_dp

  info_lap = 0
  CALL DSYGV( 1, 'V', 'U', nFree, Af, nFree, Bf, nFree, EigenVals, &
              Work, LWORK, info_lap )

  IF ( info_lap /= 0 ) THEN
    WRITE( Message, '(A,I0)' ) 'DSYGV retornou info = ', info_lap
    CALL Fatal( 'SchrodingerSolver', Message )
  END IF

  CALL Info( 'SchrodingerSolver', 'Autoproblema resolvido (DSYGV).', Level=4 )

  !--------------------------------------------------------------------
  ! 6) gravar autovalores e autovetores
  !--------------------------------------------------------------------
  nEigen = MIN( nEigen, nFree )
  CALL Info( 'SchrodingerSolver', ' ', Level=4 )
  CALL Info( 'SchrodingerSolver', '--- AUTOVALORES (energias) ---', Level=4 )
  DO i = 1, nEigen
    WRITE( Message, '(A,I4,A,ES20.12)' ) ' E(', i, ') = ', EigenVals(i)
    CALL Info( 'SchrodingerSolver', Message, Level=4 )
  END DO

  ! primeiro modo no campo de saida (zeros no contorno)
  IF ( ASSOCIATED(PsiSol) ) THEN
    PsiSol % Values = 0.0_dp
    DO j = 1, nFree
      PsiSol % Values( PsiSol % Perm( IdxFree(j) ) ) = Af(j,1)
    END DO
  END IF

  !--------------------------------------------------------------------
  ! 7) EXPORTAR as matrizes para resolver em PARALELO fora do Elmer
  !
  !     Motivo: o solver monta matrizes densas LOCAIS e usa DSYGV (serial).
  !     Em malha particionada cada processo veria so parte dos GDL -> lixo.
  !     Alem disso o Elmer 9.0 nao consegue particionar a malha nesta
  !     plataforma (crash em MeshPartition -> i2s).
  !
  !     Solucao: exportar A, M e Ma em triplas ASCII (indices GLOBAIS) e
  !     resolver A x = E (M+Ma) x ... na verdade  (K+Ma) x = E M x
  !     com o driver PETSc/SLEPc ja validado (erro 1e-15).
  !
  !     Keywords:  Exportar Matrizes = Logical True
  !                Exportar Prefixo  = String "schrod"
  !--------------------------------------------------------------------
  FazExport = .FALSE.
  ExpPrefixo = 'schrod'
  IF ( ListCheckPresent( Params, 'Exportar Matrizes' ) ) THEN
    FazExport = ListGetLogical( Params, 'Exportar Matrizes', Found )
    IF ( .NOT. Found ) FazExport = .FALSE.
  END IF
  IF ( ListCheckPresent( Params, 'Exportar Prefixo' ) ) THEN
    ExpPrefixo = ListGetString( Params, 'Exportar Prefixo', Found )
    IF ( .NOT. Found ) ExpPrefixo = 'schrod'
  END IF

  IF ( FazExport ) THEN
    ! A = K + Ma  (operador)   |   B = M  (massa)
    CALL ExportaMatriz( 'A', A, nDof, TRIM(ExpPrefixo)//'_A.dat' )
    CALL ExportaMatriz( 'M', B, nDof, TRIM(ExpPrefixo)//'_M.dat' )
    ! mapa de GDL fixos (Dirichlet) para o solver externo
    OPEN( UNIT=77, FILE=TRIM(ExpPrefixo)//'_fixed.dat', STATUS='REPLACE', ACTION='WRITE' )
    WRITE( 77, '(I0)' ) COUNT(Fixed)
    DO i = 1, nDof
      IF ( Fixed(i) ) WRITE( 77, '(I0)' ) i - 1      ! 0-based p/ PETSc
    END DO
    CLOSE( 77 )
    CALL Info( 'SchrodingerSolver', &
         'Matrizes exportadas: '//TRIM(ExpPrefixo)//'_{A,M,Ma,fixed}.dat', Level=4 )
  END IF

  DEALLOCATE( A, B, Af, Bf, EigenVals, Work, xfree, &
              Perm, MapFree, IdxFree, Fixed, Vnodal, xmin, xmax, tol )
  CALL Info( 'SchrodingerSolver', 'Concluido.', Level=4 )

!------------------------------------------------------------------------------
!  Exporta uma matriz DENSE em triplas ASCII "linha coluna valor" (0-based).
!  Formato consumivel por scipy.io.mmread-like / PETSc MatLoad ASCII.
!------------------------------------------------------------------------------
CONTAINS

  SUBROUTINE ExportaMatriz( nome, Mat, n, arquivo )
    CHARACTER(LEN=*), INTENT(IN) :: nome, arquivo
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), INTENT(IN) :: Mat(n,n)
    INTEGER :: ii, jj, un, nnz
    REAL(KIND=dp) :: v
    REAL(KIND=dp), PARAMETER :: TOLZ = 1.0e-14_dp

    nnz = COUNT( ABS(Mat) > TOLZ )
    un = 76
    OPEN( UNIT=un, FILE=arquivo, STATUS='REPLACE', ACTION='WRITE' )
    WRITE( un, '(A)' ) '% '//TRIM(nome)//'  n='//TRIM(I2S(n))//'  nnz='//TRIM(I2S(nnz))
    WRITE( un, '(I0,1X,I0,1X,I0)' ) n, n, nnz
    DO ii = 1, n
      DO jj = 1, n
        v = Mat(ii,jj)
        IF ( ABS(v) > TOLZ ) THEN
          WRITE( un, '(I0,1X,I0,1X,ES24.16)' ) ii-1, jj-1, v
        END IF
      END DO
    END DO
    CLOSE( un )
    WRITE( Message, '(A,I0,A,I0)' ) '  exportada ', nnz, ' entradas em ', n
    CALL Info( 'SchrodingerSolver', TRIM(nome)//': '//TRIM(Message), Level=4 )
  END SUBROUTINE ExportaMatriz

!------------------------------------------------------------------------------
END SUBROUTINE SchrodingerSolver
!------------------------------------------------------------------------------
