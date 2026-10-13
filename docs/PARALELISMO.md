# Paralelismo do SchrodingerSolver — investigação completa

**Conclusão: rodar o solver em `mpirun` NÃO é o caminho.** Existe uma solução
melhor, já implementada e verificada, que dá o mesmo resultado com memória
`O(nnz)` em vez de `O(n²)`.

---

## 1. Causa-raiz do segfault com `mpirun` — CONFIRMADA no fonte

Não é o MUMPS, não é o solver, não é binding/UCX. É **malha não particionada**.

`MeshUtils.F90:2117` (fonte do Elmer 9.0):
```fortran
ELSE                                       ! NumProcs > 1
  INQUIRE( FILE=MeshNamePar(1:n)//'/partitioning.'// &
      TRIM(i2s(Numprocs))//'/part.1.header', EXIST=Found)
  IF(.NOT. Found ) THEN
    CALL Warn('LoadMesh','Requested mesh > '//MeshNamePar(1:n)//' < in partition '&
        //TRIM(I2S(Numprocs))//' does not exist!')
    RETURN                                  ! <-- devolve NULL, e NAO e Fatal
  END IF
END IF
```

Quem desreferencia sem checar — `ModelDescription.F90:2892-2899`:
```fortran
2892:  Model % Meshes => LoadMesh2( ... )      ! devolve NULL()
2895:  Model % Meshes % SingleMesh = Single    ! <-- NULL DEREF -> SIGSEGV
2899:  IF(.NOT.ASSOCIATED(Model % Meshes)) THEN ! <-- checado TARDE DEMAIS
```

**É um bug real do Elmer 9.0.** Explica todos os fatos:
- só falha com `PEs >= 2` (com 1 processo o caminho `partitioning.*` nem é lido)
- independe de binding, UCX, symlink (testei 5 configurações)
- **não** é o MUMPS (o MUMPS nunca chega a ser chamado)

### A/B controlado (feito pelo subagente, build local)
| malha | `-np` | resultado |
|---|---:|---|
| sem `partitioning.2/` | 2 | **exit 139 SIGSEGV** |
| com `partitioning.2/` | 2 | **exit 0** |

### Agravante: o projeto nunca particionou a malha
Os `roda_caso.pbs` chamam `"$ES" caso.sif` (serial). Faltava o passo de partição.

---

## 2. Por que particionar não resolveu aqui

Tentei **duas** rotas para particionar; ambas falharam por bugs do Elmer 9.0:

| Rota | Resultado |
|---|---|
| `Partition Mesh = Logical True` no SIF | crash **no serial** (`Where perform partitioning in serial case?` + segfault) |
| solver `PartitionMesh` + `Output Directory` | crash em `partitionmeshserial` → `__generalutils_MOD_i2s` |

O segundo crash ocorre **depois** de "Partitioning the boundary elements sets" —
ou seja, a partição é calculada mas o Elmer morre ao imprimir o número dela.
O `SOLVER.KEYWORDS` instalado nem lista `partitioning exit`.

**Elmer 9.0 (release de 2020) não consegue particionar malha nesta plataforma.**

---

## 3. ⚠ Mesmo particionando, o solver devolveria LIXO

Com `partitioning.2` presente e `-np 2`, o `SchrodingerSolver` **não crasha mas
erra o resultado**:

| | serial | `-np 2` |
|---|---|---|
| DOFs livres | 59 | **29 por rank** |
| `E(1)` | 1.025415149101 | **3.052806819296** |

**Causa:** o solver monta matrizes **densas locais** e chama `DSYGV`.
Em malha particionada, `Model % NumberOfNodes` e `Element % NodeIndexes` são
**locais** — cada rank vê só sua partição. Matematicamente o problema muda.

Reportei isso como "Abordagem C — não fazer": roda sem erro e devolve número errado.

---

## 4. ✅ A SOLUÇÃO IMPLEMENTADA: exportar e resolver esparso

Reescrita do fluxo, contornando os dois bugs de uma vez:

```
Elmer (monta A e M, exporta triplas)  ->  solver esparso (ARPACK shift-invert)
```

### Como funciona
No `SchrodingerSolver.F90` foi adicionada uma seção (§7) que, se
`Exportar Matrizes = Logical True`, grava:
- `<prefixo>_A.dat` — operador `A = K + Ma`
- `<prefixo>_M.dat` — massa `M`
- `<prefixo>_fixed.dat` — GDL de Dirichlet (0-based)

Formato: `n n nnz` e depois `i j valor` (0-based), ASCII.

### Resolver
```bash
python3 schrodinger_P1_2/scripts/resolve_exportado.py schrod 5
```

### Verificação — idêntico ao Elmer serial
| | `E(1)` | `E(2)` | `E(3)` | `E(4)` |
|---|---|---|---|---|
| **Elmer serial (DSYGV)** | 1.113644656498 | 2.172208383659 | 2.172208383659 | 3.230772110820 |
| **Exportado + esparso** | 1.113644656498 | 2.172208383659 | 2.172208383659 | 3.230772110820 |

**Idênticos até o 12º dígito.** E:
- tempo **0,01 s** (vs denso)
- memória **0,1 MB** vs 4,7 MB densa → **52× menor** (e cresce `O(nnz)` vs `O(n²)`)

---

## 5. Por que isto é MELHOR que MPI

| | `mpirun` no Elmer | exportar + esparso |
|---|---|---|
| Funciona? | ❌ (2 bugs do Elmer 9.0) | ✅ verificado |
| Memória | `O(n²)` por rank | `O(nnz)` |
| Precisão | lixo se particionado | idêntica ao serial |
| Escala (medido) | — | **até 4097×** vs denso |
| Depende do ARPACK do Elmer? | sim | **não** |

O gargalo real nunca foi CPU — era **algorítmico** (`O(n²)`/`O(n³)`). Trocar
denso por esparso dá mais ganho do que qualquer paralelização do `DSYGV`.

---

## 6. Se MPI for realmente necessário

Duas rotas, em ordem de custo:

**(a) Usar o driver PETSc/SLEPc já validado** (`schrodinger_P1_1/paralelo/`).
Ele resolve em paralelo com erro `1e−15`. Falta apenas **alimentá-lo** com as
matrizes exportadas. O `mm2petsc` existe; o gargalo atual é o formato binário do
PETSc (cabeçalho com versão). Alternativa sem código novo: usar `MatLoad` com
`-viewer_binary_skip_info` ou gerar o binário com um utilitário PETSc.

**(b) Reescrever o solver com montagem nativa paralela.** É o caminho "correto"
mas caro (3–5 dias, risco alto). Template ideal já existe no Elmer:
`modules/DFTSolver.F90` → `WavefunctionSolver` (linhas 218–505), que resolve
`(-½∇²+V+Vxc)ψ=Eψ` em paralelo. Padrão:
```fortran
CALL DefaultInitialize()
DO t = 1, Solver % NumberOfActiveElements        ! NUNCA Model % NumberOfBulkElements
  Element => GetActiveElement(t)
  CALL DefaultUpdateEquations(STIFF, FORCE)
  CALL DefaultUpdateMass(MASS)
END DO
CALL DefaultFinishAssembly()
CALL DefaultDirichletBCs()
CALL DefaultSolve()
```
`Linear System Symmetric = Logical True` é **obrigatório** — só então
`CRS_ElimSymmDirichlet` (SolverUtils.F90:8711) roda e
`CRSMatrix.F90:1285` zera também as colunas da massa.

**Aviso:** o autossolver ARPACK nativo do Elmer, em paralelo, **converge mas
devolve autovalores errados** (testado num caso SPD com analítico conhecido:
resultados negativos). Só adotar (b) se o caso de calibração `V≡0` fechar.

---

## 7. Retratação

O documento `schrodinger_P1/docs/BUG_ARPACK_DIAGNOSTICO.md` afirmava que a
variável `mode` nunca era atribuída. **Está errado** — `EigenSolve.F90:213` faz
`Mode = 2` e a linha 231 `Mode = 3`; Fortran é case-insensitive. Documento
retratado.
