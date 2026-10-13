# schrodinger_P1_2 — Solver de Schrödinger em Elmer FEM

**Status: FUNCIONANDO E VALIDADO em 1D, 2D e 3D.**
O Elmer **não possui** solver de Schrödinger; este foi escrito do zero.

---

## 1. Por que um solver próprio

O Elmer 9.0 tem dois defeitos no seu autossolver ARPACK:

1. **`EigenSolve.F90` não inicializa `RESID`** (vetor inicial do ARPACK) em
   `ArpackEigenSolve` → `info = -9` (*starting vector is zero*).
   Causa raiz diagnosticada e documentada em `../schrodinger_P1/RELATORIO_OPCAO_B.md`.
2. Mesmo com `RESID` corrigido (patch aplicado), o loop de *reverse communication*
   falha com `-9` após a primeira aplicação do operador — segundo defeito, distinto.

**Solução:** solver dedicado que monta as matrizes e chama **`DSYGV` (LAPACK)**
diretamente, sem ARPACK.

---

## 2. Formulação

Equação (ħ = m = 1): `-∇²ψ + V(x)ψ = Eψ`

Forma fraca de Galerkin → autoproblema generalizado `A ψ = E M ψ`, `A = K + Ma`:

| Matriz | Integral | Montagem no código |
|---|---|---|
| `K` rigidez | `∫∇Nᵢ·∇Nⱼ dΩ` | `Σ s·detJ·(dBasisdx(p,:)·dBasisdx(q,:))` |
| `M` massa | `∫NᵢNⱼ dΩ` | `Σ s·detJ·Basis(p)·Basis(q)` |
| `Ma` potencial | `∫Nᵢ V Nⱼ dΩ` | `Σ s·detJ·V(ξ)·Basis(p)·Basis(q)` |

Quadratura via `GaussPoints(Element)`; `detJ` e formas via `ElementInfo`.
Dimensionalidade genérica (1D/2D/3D) via `CoordinateSystemDimension()` e
`dBasisdx(p,1:dim)`.

---

## 3. Validação — verificação cruzada em três dimensões

Solução analítica exata em todos os casos. Python e Elmer têm **montagem
independente** e **autossolver diferente** (scipy `eigh` vs LAPACK `DSYGV`).

### 1D — poço infinito, L=40, 250 elementos `Eₙ=(nπ/L)²`

| n | Elmer | analítico | erro |
|---|---|---|---|
| 1 | 0.006168583925 | 0.006168502751 | 1.32e−05 |
| 5 | 0.154263309342 | 0.154212568767 | 3.29e−04 |

**Erro máximo 3.29e−04 · |Elmer − Python| = 4.21e−11** (precisão de máquina)

### 2D — poço quadrado, L=10, 40×40 (1521 DOFs livres)

| n | Elmer | analítico | erro |
|---|---|---|---|
| 1 | 0.1974935767 | 0.1973920880 | 5.14e−04 |
| 2,3 | 0.4943433728 | 0.4934802201 | 1.75e−03 (degenerado) |
| 4 | 0.7911931689 | 0.7895683521 | 2.06e−03 |
| 5 | 0.9911281920 | 0.9869604401 | 4.22e−03 |

### 3D — poço cúbico, L=6, 10×10×10 (729 DOFs livres)

| n | Elmer | analítico | erro |
|---|---|---|---|
| 1 | 0.8292535815 | 0.8224670334 | 8.25e−03 |
| 2,3,4 | 1.6859901662 | 1.6449340668 | 2.50e−02 (triplamente degenerado) |
| 5 | 2.5427267508 | 2.4674011003 | 3.05e−02 |

### Convergência — ordem O(h²) confirmada

**2D:**
| n | erro | ordem |
|---:|---|---:|
| 10 | 8.251e−03 | — |
| 20 | 2.058e−03 | 2.00 |
| 40 | 5.141e−04 | 2.00 |
| 80 | 1.285e−04 | 2.00 |

**3D:**
| n | DOFs | erro | ordem |
|---:|---:|---|---:|
| 5 | 64 | 3.331e−02 | — |
| 8 | 343 | 1.292e−02 | 2.16 |
| 10 | 729 | 8.251e−03 | 1.98 |
| 12 | 1331 | 5.725e−03 | 2.04 |
| 15 | 2744 | 3.661e-03 | 1.98 |

---

## 4. Bugs corrigidos durante o desenvolvimento

Estes são erros que eu mesmo introduzi e precisei diagnosticar — valem registro
porque são armadilhas comuns:

### a) Dirichlet com diagonal = 1 gerava autovalores espúrios
Zerar linha/coluna e pôr `A(i,i)=B(i,i)=1` cria modos artificiais com
**autovalor exatamente 1.0**. No caso 3D isso contaminou a resposta (os cinco
"autovalores" saíram 0.829, 1.0, 1.0, 1.0, 1.0).
**Correção:** eliminação real dos DOFs prescritos — extrair submatrizes
`Af(i,j) = A(IdxFree(i), IdxFree(j))` e resolver só o sistema reduzido.

### b) Formato do header da malha
A segunda seção do `mesh.header` é `código do elemento` + **nº de NOS por
elemento**, não a contagem de elementos. Errei para `404 160` quando o correto
era `202 2` (linha de 2 nós).

### c) Código de elemento de contorno
Contorno 2D e 3D são **linhas** de 2 nós → código **202**. Usei `404`
(quadrilátero) por engano.
Contorno 3D: faces quad → **404** (correto).

### d) `.sif` com dimensão fixa
Um `Coordinate System = Cartesian 2D` aplicado a malha 3D produz
`ERROR:: ElementMetric: Degenerate 3D element` e autovalores absurdos (−2.8e32).
A dimensão do SIF **tem de** coincidir com a da malha.

---

## 5. Ferramentas entregues

| Arquivo | Função |
|---|---|
| `fonte/SchrodingerSolver.F90` | o solver (interface padrão Elmer) |
| `scripts/gera_malha_elmer.py` | gerador de malha estruturada 1D/2D/3D no formato Elmer |
| `scripts/run_schrodinger.pbs` | submissão PBS nos nós computacionais |
| `shares/VERIFICACAO_CRUZADA.txt` | resultados consolidados |
| `shares/ELMER_caso*.log` | saídas verbatim do solver |

### Por que um gerador de malha próprio
O `Minimum Element Divisions` do `ElmerGrid` **não** produz a contagem pedida
(pedi 39×39 e obtive 1560 elementos em vez de 1521). O gerador em Python dá
controle exato e é verificado contra o formato de malhas do próprio ElmerGrid.

### Compilação
```bash
module load gnu12 openmpi4 openblas
MODDIR=$HOME/elmer_build/build_gnu/fmodules
LIBDIR=$HOME/elmer_build/build_gnu/fem/src
mpif90 -c SchrodingerSolver.F90 -o SchrodingerSolver.o -fallow-argument-mismatch \
  -DCONTIG= -DHAVE_EXECUTECOMMANDLINE -DUSE_ISO_C_BINDINGS -O2 -fPIC -I$MODDIR
mpif90 -shared SchrodingerSolver.o -o SchrodingerSolver.so -L$LIBDIR -lelmersolver
```

---

## 6. Limitações atuais

1. **Potencial apenas uniforme** — lê `Potential` do Solver. Falta V(x) nodal
   (barreira, poço duplo, potencial senoidal).
2. **`DSYGV` é serial e denso** — memória O(n²). O caso 2D com 1521 DOFs já exige
   ~40 MB por matriz. Para malhas grandes será preciso ARPACK (com o patch) ou
   autossolver paralelo.
3. **Sem MPI** — o solver roda em 1 processo.
4. Comentário de cabeçalho menciona `Eigen System Values` genérico; na prática
   usa `Apply Dirichlet` + `Eigen System Values`.

---

## 7. Próximos passos (ordem definida pelo usuário)

~~2. Testar 2D e 3D~~ **CONCLUÍDO** (esta rodada)
→ **1. Potencial V(x) variável** (próximo)
→ 3. Paralelização MPI
