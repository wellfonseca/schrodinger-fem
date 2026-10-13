# Tutorial — Poço duplo — tunelamento quântico

**Directory:** `1D_poco_duplo`
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
| Potencial | V = 0.01(x² − 4)² |

O problema a resolver é a equação de Schrödinger estacionária

    −∇²ψ + V(r) ψ = E ψ

cuja forma fraca leva ao autoproblema generalizado `A ψ = E M ψ`, com
`A = K + Mᵃ`. Os detalhes da formulação estão no
`Manual_Modelo_Schrodinger.md`.

### Solução analítica

    splitting de tunelamento: ΔE = E₁ − E₀ > 0

**Valor de referência:** (sem forma fechada simples) — verifica-se ΔE > 0 e simetria par/ímpar

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
    Real MATC "0.01*(tx*tx-4.0)*(tx*tx-4.0)"
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
| 1 | `2.241105142386E−01` | (par) | — |
| 2 | `6.000280051832E−01` | (ímpar) | — |
| 3 | `1.274687373480E+00` | — | — |
| 4 | `2.059803060096E+00` | — | — |

**Esperado:** E₀ ≈ 0.2241, E₁ ≈ 0.6000  →  ΔE ≈ 0.3758

### Verificação do potencial

O solver reporta no log os extremos do campo `V` lido:

```
SchrodingerSolver: V(x) lido como campo nodal de Body Force.
SchrodingerSolver:   V_min = ...
SchrodingerSolver:   V_max = ...
```

Confira esses valores **antes** de confiar nos autovalores.

### Observações

Potencial com **dois mínimos** em `x = ±2`, separados por uma barreira em
`x = 0` com altura `V(0) = 0.01·16 = 0.16`.

O fenômeno físico é o **tunelamento**: classicamente uma partícula com
`E < 0.16` ficaria presa em um dos poços, mas quanticamente a função de onda
se espalha pelos dois. O resultado é que o nível fundamental se **divide em
dois**, com splitting `ΔE = E₁ − E₀ > 0`.

**Como verificar:** o splitting deve ser **positivo** e **diminuir** se a
barreira for aumentada (aumentar o coeficiente de `x²` ou reduzir `x⁴`).

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

Este é o caso mais interessante visualmente:

- **modo 1** (E₀): função **simétrica** — duas "lombadas" de mesmo sinal, uma em
  cada poço. Estado ligante (bonding).
- **modo 2** (E₁): função **antissimétrica** — as duas lombadas com sinais
  **opostos**. Estado antiligante (antibonding).

A diferença `E₁ − E₀` é a **energia de splitting**, diretamente ligada à
probabilidade de tunelamento. No ParaView, selecione `psi` e observe a
alternância de sinal entre as duas cavidades.

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
