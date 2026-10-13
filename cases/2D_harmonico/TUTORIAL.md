# Tutorial — Oscilador harmônico isotrópico em 2D

**Directory:** `2D_harmonico`
**Solvers:** SchrodingerSolver, ResultOutputSolve
**Tools:** ElmerSolver, ParaView
**Dimensions:** 2D, Estacionário
**Tipo:** autoproblema de valores próprios (generalized eigenvalue problem)

---

## 1. Definição do caso

O domínio é `[−L/2, +L/2]²` com

| Parâmetro | Valor |
|---|---|
| Comprimento `L` | `20.0` |
| Elementos por direção | `40` |
| Origem | `-10.0` |
| Potencial | V = 0.25(x² + y²) |

O problema a resolver é a equação de Schrödinger estacionária

    −∇²ψ + V(r) ψ = E ψ

cuja forma fraca leva ao autoproblema generalizado `A ψ = E M ψ`, com
`A = K + Mᵃ`. Os detalhes da formulação estão no
`Manual_Modelo_Schrodinger.md`.

### Solução analítica

    E_n = 2√α (n + d/2),   α = 0.25, d = 2  →  E_n = n + 1

**Valor de referência:** E₀ = 1.0; E₁ = E₂ = 2.0 (duplo); E₃ = E₄ = 3.0 (duplo)

---

## 2. Procedimento de solução

> ⚠ **Toda simulação deve rodar nos nós computacionais, nunca no headnode.**

### 2.1 Gerar a malha

A malha é gerada pelo script `gerar_malha.py`, que escreve o formato nativo do
Elmer (`mesh.header`, `mesh.nodes`, `mesh.elements`, `mesh.boundary`):

```bash
python3 gerar_malha.py
```

### 2.2 O arquivo de entrada (SIF)

**`Header`** e **`Simulation`**:

```
Header
  Mesh DB "." "malha"
End

Simulation
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Max Output Level = 5
End
```

O `Coordinate System` **deve** corresponder à dimensão da malha.

**`Body`, `Material` e `Body Force`**:

```
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End

Material 1
  Density = 1.0
End

Body Force 1
  Potential = Variable Coordinate 1, Coordinate 2
    Real MATC "0.25*(tx*tx + tx(1)*tx(1))"
End
```

**`Equation`** — dois solvers ativos:

```
Equation 1
  Active Solvers(2) = 1 2
End
```

**`Solver 1`** — o solver de Schrödinger:

```
Solver 1
  Equation = "Schrodinger"
  Procedure = "SchrodingerSolver" "SchrodingerSolver"
  Variable = "psi"
  Variable DOFs = 1
  Eigen System Values = 5
  Apply Dirichlet = Logical True
  Linear System Solver = Direct
  Linear System Direct Method = "umfpack"
End
```

`Eigen System Values = 5` pede os **cinco menores** autovalores.
`Apply Dirichlet = Logical True` impõe `ψ = 0` no contorno, detectado
automaticamente pela geometria — **não** é preciso declarar `Boundary Condition`.

**`Solver 2`** — a saída para ParaView:

```
Solver 2
  Exec Solver = "After Simulation"
  Equation = "ResultOutput"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = File "run"
  Output Format = String "vtu"
End
```

### 2.3 Executar

Via PBS (recomendado — roda num nó computacional):

```bash
qsub roda_caso.pbs
```

Interativamente, dentro de um nó (`qsub -I`):

```bash
$ES caso.sif > run.log 2>&1
```

---

## 3. Resultados

| Modo | `E` calculado | Referência | Erro relativo |
|---|---|---|---|
| 1 | `1.028585384356E+00` | 1.0 | 2.86e−02 |
| 2 | `2.043966018888E+00` | 2.0 | 2.20e−02 |
| 3 | `2.043966018888E+00` | 2.0 | 2.20e−02 |
| 4 | `3.059346653419E+00` | 3.0 | 1.98e−02 |

**Esperado:** E₀ ≈ 1.0286; degenerescência DUPLA em E₂ = E₃

### Verificação do potencial

O solver reporta no log os extremos do campo `V` lido:

```
SchrodingerSolver: V(x) lido como campo nodal de Body Force.
SchrodingerSolver:   V_min = ...
SchrodingerSolver:   V_max = ...
```

Confira esses valores **antes** de confiar nos autovalores.

### Observações

**Este caso introduz a sintaxe do MATC multidimensional.** O separador de
variáveis é **vírgula**:

    Potential = Variable Coordinate 1, Coordinate 2

e as componentes se acessam como `tx` (primeira) e `tx(1)` (segunda).
Usar **espaço** (`Coordinate 1 Coordinate 2`) faz o Elmer interpretar tudo como
**um único nome** inexistente.

**Verificação do potencial:** `V_min = 0`, `V_max = 0.25·(10²+10²) = 50.0`.

**Degenerescência dupla em E₂ = E₃:** os estados `(nₓ,n_y) = (1,2)` e `(2,1)`
têm a mesma energia por simetria da permutação. Isso só aparece se a montagem e
a simetria do domínio estiverem corretas — é um teste independente e forte.

---

## 4. Visualização

O caso gera o arquivo VTU (malha + o campo `psi`, autofunção do modo 1):

```bash
paraview malha/run_t0001.vtu
```

> O Elmer escreve o `.vtu` no **diretório da malha**, por isso o caminho é
> `malha/run_t0001.vtu` e não `run.vtu` na raiz.

No ParaView: selecione `psi` em **Coloring** para ver a autofunção. Use os
filtros *Contour* ou *Slice* (em 3D) para inspecionar a estrutura nodal.

### O que observar

O estado fundamental é uma **gaussiana 2D** centrada na origem, com simetria
radial. No ParaView, o modo 1 deve parecer um "morro" circular e suave.

O primeiro estado excitado é **duplamente degenerado**: existem duas
autofunções independentes com a mesma energia, uma com um nó vertical e outra
com um nó horizontal (análogas a `x·e^{−r²/2}` e `y·e^{−r²/2}`). Como o solver
devolve apenas uma combinação (o modo 1 no `.vtu`), para ver as duas seria
necessário resolver as matrizes exportadas.

---

## 5. Comandos completos

O arquivo `command.txt`, nesta mesma pasta, tem todos os comandos de terminal
passo a passo (do `ssh` ao `paraview`), incluindo as variáveis de ambiente
necessárias.

---

## 6. Referências

| Documento | Conteúdo |
|---|---|
| `Manual_Modelo_Schrodinger.md` | física, keywords, convenções de unidades |
| `Manual_Solver_Schrodinger.md` | implementação, compilação, diagnóstico |
| `docs/VALIDACAO_POTENCIAL.md` | validação do potencial `V(x)` |
| `docs/VALIDACAO_V_2D3D.md` | validação 2D/3D e convergência `O(h²)` |
| `docs/SINTAXE_MATC_POTENCIAL.md` | sintaxe do MATC (separador vírgula) |
| `docs/PARALELISMO.md` | por que não usar `mpirun` com este solver |
| `../LEIA-ME.md` | visão geral dos casos |
| `../../RESULTADOS_TODOS_CASOS.md` | resultados de todos os 12 casos |

---

*Última atualização: 08/10/2026*
