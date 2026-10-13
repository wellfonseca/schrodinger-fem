# Tutorial — Oscilador harmônico em 1D

**Directory:** `1D_harmonico`
**Solvers:** SchrodingerSolver, ResultOutputSolve
**Tools:** ElmerSolver, ParaView
**Dimensions:** 1D, Estacionário
**Tipo:** autoproblema de valores próprios (generalized eigenvalue problem)

---

## 1. Definição do caso

O domínio é `[−L/2, +L/2]` com

| Parâmetro | Valor |
|---|---|
| Comprimento `L` | `20.0` |
| Elementos por direção | `400` |
| Origem | `-10.0` |
| Potencial | V = x²  (oscilador harmônico isotrópico) |

O problema a resolver é a equação de Schrödinger estacionária

    −∇²ψ + V(r) ψ = E ψ

cuja forma fraca leva ao autoproblema generalizado `A ψ = E M ψ`, com
`A = K + Mᵃ`. Os detalhes da formulação estão no
`Manual_Modelo_Schrodinger.md`.

### Solução analítica

    E_n = 2√α · (n + ½),   com α = 1

**Valor de referência:** E₀ = 1.0,  E₁ = 3.0,  E₂ = 5.0  (espaçamento 2.0)

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
  Coordinate System = Cartesian 1D
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
  Potential = Variable Coordinate 1
    Real MATC "1.0*tx*tx"
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
| 1 | `1.000572892258E+00` | 1.0 | 5.73e−04 |
| 2 | `3.001197696943E+00` | 3.0 | 3.99e−04 |
| 3 | `5.002447062087E+00` | 5.0 | 4.89e−04 |
| 4 | `7.004320694552E+00` | 7.0 | 6.17e−04 |
| 5 | `9.006818301114E+00` | 9.0 | 7.57e−04 |

**Esperado:** E₀ ≈ 1.0006, E₁ ≈ 3.0012, espaçamento ≈ 2.0

### Verificação do potencial

O solver reporta no log os extremos do campo `V` lido:

```
SchrodingerSolver: V(x) lido como campo nodal de Body Force.
SchrodingerSolver:   V_min = ...
SchrodingerSolver:   V_max = ...
```

Confira esses valores **antes** de confiar nos autovalores.

### Observações

**Este é o caso que determina a convenção de unidades do solver.** O poço
infinito com `L = 40` **não** distingue as convenções `−ψ''` e `−½ψ''`; o
oscilador harmônico distingue, porque o autovalor depende de `√α` e não de
`1/L²`.

Com `V = x²` e a convenção implementada (`−ψ'' + Vψ = Eψ`), o resultado é
`E_n = 2(n + ½)`. A forma usual `−½ψ'' + Vψ = Eψ` daria `E_n = (n + ½)` —
**metade**. Ver `Manual_Modelo_Schrodinger.md` §6.

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

A autofunção do estado fundamental é uma **gaussiana** centrada em `x = 0`:

    ψ₀(x) ∝ exp(−x²/2)

No ParaView, o modo 1 deve aparecer como um pico suave e simétrico no centro.
O modo 2 deve ter **um nó** exatamente em `x = 0` (função ímpar).

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
