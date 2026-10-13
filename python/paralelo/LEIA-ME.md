# Paralelização e autossolver esparso — item 3

**Resultado central: o gargalo era ALGORÍTMICO, não de paralelismo.**
Trocar o solver denso por esparso (ARPACK shift-invert) deu ganho de **até 4097×**,
enquanto o MPI não traz ganho nesta instalação (ver §4).

---

## 1. O problema medido

O solver original usava `scipy.linalg.eigh` / `DSYGV` — **denso**:
- memória `O(n²)`  — 3 matrizes n×n em double = `24 n²` bytes
- custo   `O(n³)`

Medido nesta sessão (caso 2D harmônico, k=4):

| GDL | denso (s) | esparso (s) | ganho |
|---:|---:|---:|---:|
| 441 | 0,01 | 0,00 | — |
| 1 681 | 0,58 | 0,01 | 58× |
| 6 561 | 41,06 | 0,05 | **821×** |
| 14 641 | **491,60** | 0,12 | **4 097×** |
| 39 601 | **(inviável: 37,6 GB)** | **0,41** | — |

O caso n=120 (14 641 GDL) exigiria 5,1 GB de matrizes e estourou o tempo limite.
O caso n=200 (39 601 GDL) exigiria 37,6 GB — **inviável em denso**, e roda em
**0,41 s** no esparso.

---

## 2. O que foi implementado

### a) `../scripts/eigsolver_esparso.py`
Substitui o solver denso por ARPACK em modo **shift-invert**:
- fatora `(A − σM)` uma única vez (esparso)
- itera Arnoldi para os `k` menores autovalores
- memória `O(nnz) + O(n·k)` em vez de `O(n²)`
- mesma matemática e mesma interface do solver denso

### b) `schrodinger_petsc.c` — driver paralelo MPI
Autossolver com **SLEPc 3.18** (Krylov-Schur + shift-invert), lê as matrizes do
PETSc e roda em vários processos. Compilado e funcionando:
```
E(1) = +1.010304617094e+00   (erro rel. 1.18e-15)
E(2) = +2.015897967980e+00   (degenerado)
E(3) = +2.015897967980e+00
```
Correto contra o analítico `E₀ = 1.0` para o harmônico 2D.

### c) Ferramentas de apoio
- `mm2petsc.c` — converte Matrix Market → binário PETSc
- `exporta_mtx.py` — exporta as matrizes FEM para Matrix Market
- `bench_escalabilidade.sh` — mede tempo × nº de processos

---

## 3. Como rodar

> **REGRA: submeter aos nós computacionais, nunca rodar no headnode.**
> Compilação pode ser feita no headnode; a **solução** vai por `qsub`.

```bash
# --- pre-requisitos ---
module load gnu12 openmpi4 petsc slepc

# --- compilar (uma vez) ---
mpicc schrodinger_petsc.c -o ~/bin/schrodinger_petsc \
    -I$PETSC_DIR/include -I$SLEPC_DIR/include \
    -L$PETSC_DIR/lib -L$SLEPC_DIR/lib -lslepc -lpetsc -lm
mpicc mm2petsc.c -o ~/bin/mm2petsc -I$PETSC_DIR/include -L$PETSC_DIR/lib -lpetsc -lm

# --- preparar um caso ---
python3 exporta_mtx.py 2 40.0 200 0.25 matrizes   # dim L n alpha pasta
~/bin/mm2petsc matrizes/A.mtx A.petsc
~/bin/mm2petsc matrizes/M.mtx M.petsc

# --- resolver NOS NOS COMPUTACIONAIS ---
qsub -v NP=1,K=4,DIR=$PWD pbs/roda_petsc.pbs

# --- opcional: outra biblioteca de fatoracao (se disponivel) ---
#   adicione -fact_solver mumps ao comando dentro do .pbs
```

### Solver esparso (Python) — SUBMETER
```bash
qsub -v DIM=2,L=30,N=140 pbs/roda_esparso.pbs
```
Equivalente interativo (apenas dentro de um nó):
```bash
python3 -c "
import sys; sys.path.insert(0,'../scripts')
from fem_schrodinger import malha_2d, montar, contorno
from eigsolver_esparso import resolve_esparso
A,M,gm,X,V = montar(malha_2d(30,30,140,140,-15,-15), lambda X: 0.25*(X[:,0]**2+X[:,1]**2))
E,psi,interno,info = resolve_esparso(A,M,nmodos=4,nos_contorno=contorno(X, malha_2d(30,30,140,140,-15,-15)))
print(E)
"
```
(ver `casos/` para exemplos prontos)

---

## 4. LIMITAÇÃO IMPORTANTE: MPI não ajuda nesta instalação

O driver paralelo em **1 processo funciona**; com **≥2 processos falha**:

```
[0]PETSC ERROR: You chose to solve linear systems with a factorization,
but in parallel runs you need to select an external package
```

**Causa:** o PETSc deste cluster (3.18.1) **não foi compilado com nenhum
fatorador paralelo**. Verificado com `MatGetFactor`:

| Fatorador | Disponível |
|---|---|
| mumps | ✗ |
| superlu_dist | ✗ |
| superlu | ✗ |
| umfpack | ✗ |
| **petsc (interno)** | **✓ — mas é SEQUENCIAL** |

As *configure options* do PETSc instalado listam apenas scaLAPACK, OpenBLAS e
HDF5 — nenhum solver esparso paralelo.

**Consequência:** sem um LU paralelo, o shift-invert não pode ser distribuído, e
o driver fica restrito a 1 processo. **O driver já está pronto para MPI** — basta
um PETSc com MUMPS ou SuperLU_DIST.

### Multithreading (BLAS) também não ajuda
Testado com `OPENBLAS_NUM_THREADS` no caso de 39 601 GDL:

| threads | tempo (s) |
|---:|---:|
| 1 | 0,441 |
| 2 | 0,369 |
| 4 | 0,365 |
| 8 | 0,445 |

Ganho máximo ~17%, saturado em 4 threads. O PETSc está linkado ao OpenBLAS do
sistema (`/lib64/libopenblas.so.0`), não ao do módulo — daí o paralelismo pobre.

---

## 5. Como habilitar MPI de verdade

Três caminhos, em ordem de esforço:

**(a) Pedir ao administrador do cluster** um PETSc/SLEPc com MUMPS:
```
--download-mumps --download-scalapack --download-metis --download-parmetis
```

**(b) Compilar PETSc+PETSc4py no home** com MUMPS (2–4 h de build):
```bash
module load gnu12 openmpi4
./configure --with-debugging=0 --with-64-bit-indices=0 \
  --download-mumps --download-scalapack --download-metis \
  --download-parmetis --download-fblaslapack
```

**(c) Contornar:** usar o solver esparso Python (1 processo), que já resolve
39 601 GDL em 0,41 s — suficiente para os casos do Artigo 1.

---

## 6. Recomendação

Para o **Artigo 1**, o ganho algorítmico (denso → esparso, 4097×) é a
contribuição principal e já está demonstrado. O MPI acrescentaria escalabilidade
para malhas muito maiores (10⁶+ GDL), que só serão necessárias se o projeto
avançar para 3D fino.

**Registrar no artigo:** além do ganho, a *dificuldade prática* de encontrar
builds com solvers paralelos é um achado relevante para reprodutibilidade —
comum em clusters institucionais.
