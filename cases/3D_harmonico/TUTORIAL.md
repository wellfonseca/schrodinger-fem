# Tutorial — Oscilador harmônico isotrópico em 3D

**Directory:** `3D_harmonico`
**Solvers:** SchrodingerSolver, ResultOutputSolve
**Tools:** ElmerSolver, ParaView
**Dimensions:** 3D, Estacionário
**Tipo:** autoproblema de valores próprios (generalized eigenvalue problem)

---

## 1. Definição do caso

O domínio é `[−L/2, +L/2]³` com

| Parâmetro | Valor |
|---|---|
| Comprimento `L` | `12.0` |
| Elementos por direção | `14` |
| Origem | `-6.0` |
| Potencial | V = 0.25(x² + y² + z²) |

O problema a resolver é a equação de Schrödinger estacionária

    −∇²ψ + V(r) ψ = E ψ

cuja forma fraca leva ao autoproblema generalizado `A ψ = E M ψ`, com
`A = K + Mᵃ`. Os detalhes da formulação estão no
`Manual_Modelo_Schrodinger.md`.

### Solução analítica

    E_n = 2√α (n + d/2),   α = 0.25, d = 3  →  E_n = n + 1.5

**Valor de referência:** E₀ = 1.5; E₁ = E₂ = E₃ = 2.5 (TRIPLO)

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
  Coordinate System = Cartesian 3D
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
  Potential = Variable Coordinate 1, Coordinate 2, Coordinate 3
    Real MATC "0.25*(tx*tx + tx(1)*tx(1) + tx(2)*tx(2))"
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
| 1 | `1.625507320661E+00` | 1.5 | 8.37e−02 |
| 2 | `2.669307913349E+00` | 2.5 | 6.77e−02 |
| 3 | `2.669307913349E+00` | 2.5 | 6.77e−02 |
| 4 | `2.669307913349E+00` | 2.5 | 6.77e−02 |

**Esperado:** E₀ ≈ 1.6255; degenerescência TRIPLA em E₂ = E₃ = E₄

### Verificação do potencial

O solver reporta no log os extremos do campo `V` lido:

```
SchrodingerSolver: V(x) lido como campo nodal de Body Force.
SchrodingerSolver:   V_min = ...
SchrodingerSolver:   V_max = ...
```

Confira esses valores **antes** de confiar nos autovalores.

### Observações

Caso **3D** com potencial. Três componentes acessadas como `tx`, `tx(1)`,
`tx(2)` no MATC.

**Verificação do potencial:** `V_min = 0`,
`V_max = 0.25·(6²+6²+6²) = 27.0`.

**Degenerescência TRIPLA em E₂ = E₃ = E₄:** os estados `(1,1,2)`, `(1,2,1)` e
`(2,1,1)` têm a mesma energia. A multiplicidade **três** é exatamente o que se
espera em 3D, e distingue este caso do 2D (onde é dois). É o teste de simetria
mais forte da bateria.

**Erros maiores (6–8%):** a malha é grosseira (14³) e o estado excitado se
estende mais que o domínio. Refinar a malha ou aumentar `L` reduz o erro.

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

O estado fundamental é uma gaussiana **esfericamente simétrica**. No ParaView,
usando o filtro *Contour* ou *Slice*, o modo 1 deve mostrar simetria radial.

O primeiro estado excitado tem **degenerescência tripla** — três autofunções
independentes com a mesma energia, com nós em `x`, `y` ou `z`. A visualização
de uma única combinação mostrará apenas um dos padrões nodais.

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
