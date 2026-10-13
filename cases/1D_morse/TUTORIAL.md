# Tutorial — Potencial de Morse (potencial molecular)

**Directory:** `1D_morse`
**Solvers:** SchrodingerSolver, ResultOutputSolve
**Tools:** ElmerSolver, ParaView
**Dimensions:** 1D, Estacionário
**Tipo:** autoproblema de valores próprios (generalized eigenvalue problem)

---

## 1. Definição do caso

O domínio é `[−L/2, +L/2]` com

| Parâmetro | Valor |
|---|---|
| Comprimento `L` | `10.0` |
| Elementos por direção | `400` |
| Origem | `0.0` |
| Potencial | V = Dₑ(e^{−2β(x−rₑ)} − 2e^{−β(x−rₑ)}),  Dₑ=10, β=0.5, rₑ=3 |

O problema a resolver é a equação de Schrödinger estacionária

    −∇²ψ + V(r) ψ = E ψ

cuja forma fraca leva ao autoproblema generalizado `A ψ = E M ψ`, com
`A = K + Mᵃ`. Os detalhes da formulação estão no
`Manual_Modelo_Schrodinger.md`.

### Solução analítica

    mínimo do potencial = −Dₑ  (verificação do campo V)

**Valor de referência:** V_min = −Dₑ = −10.0 exatamente, em x = rₑ = 3

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
    Real MATC "10.0*(exp(-2*0.5*(tx-3.0)) - 2*exp(-0.5*(tx-3.0)))"
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
| 1 | `−8.481033786641E+00` | (ponto zero) | — |
| 2 | `−5.818565431911E+00` | — | — |
| 3 | `−3.656064809003E+00` | — | — |
| 4 | `−1.991029948157E+00` | — | — |

**Esperado:** V_min = −10.000000 e E₀ ≈ −8.48 (energia de ponto zero ≈ 1.52)

### Verificação do potencial

O solver reporta no log os extremos do campo `V` lido:

```
SchrodingerSolver: V(x) lido como campo nodal de Body Force.
SchrodingerSolver:   V_min = ...
SchrodingerSolver:   V_max = ...
```

Confira esses valores **antes** de confiar nos autovalores.

### Observações

O potencial de Morse é o caso de teste **assimétrico** com `V ≠ 0`: o poço não
é simétrico em torno da origem, então exercita a integração correta de `V` nos
pontos de Gauss.

**O valor de referência é o mínimo do potencial:** `V_min = −Dₑ = −10`
exatamente. O solver reporta `V_min`/`V_max` no log — é a verificação mais
direta de que o campo `V` foi lido e avaliado corretamente.

O `L = 10` com origem em `x = 0` posiciona o mínimo do poço em `x = 3`, bem
dentro do domínio, com `ψ` decaindo até zero nas duas paredes.

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

A autofunção do modo 1 fica **localizada em torno de x = 3** — o fundo do poço —
e decai exponencialmente para os dois lados. Diferente de um poço simétrico,
o pico **não** está no centro do domínio.

Este é o análogo quântico da vibração molecular: `E₀` acima do fundo é a
**energia de ponto zero**, `≈ 1.519` (comparável a `ω/2` na aproximação
harmônica local, coerente com o espaçamento observado entre níveis).

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
