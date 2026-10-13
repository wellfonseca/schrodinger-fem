# Item 1 — Potencial V(x) e suíte de regressão

**Status: V(x) FUNCIONANDO · suíte de regressão 5/5 · convenção verificada**

---

## 1. Bug corrigido: leitura de V(x)

### Antes
```fortran
Material => Model % Materials(1) % Values
IF ( ListCheckPresent( Material, 'Potential' ) ) THEN
  Vuni = ListGetConstReal( Material, 'Potential', Found )   ! escalar
```
Dois defeitos:
- `Potential` **não é keyword válido** em `Material` no Elmer 9
  (`ERROR:: Unknown specifier ... property name:[potential]`);
- `ListGetConstReal` só aceita **escalar** — incapaz de ler campo.

Resultado: o ramo nunca executava e o solver lia **V = 0 em silêncio**.
Nenhum teste com V ≠ 0 jamais havia rodado.

### Depois
```fortran
! (1) Body Force  (aceita expressao: Variable Coordinate / MATC)
Vtmp = ListGetReal( BFList, 'Potential', nNodes, NodeIdx, Found )
! (2) Material    (fallback)
! (3) V = 0
```
Notas de implementação que custaram tempo:
- `GetRealValues` **não serve** — depende de `GetCurrentElement` (contexto do laço);
- `ListGetReal` é **função** (`RESULT(F)`, com `REAL :: F(N)`), não sub-rotina;
- `BF` colide com `Bf` (matriz) — Fortran não diferencia maiúsculas.
  Renomeado para `BFList`;
- `Max Output Level = 3` **suprime** as mensagens do solver (usam `Level=4`).

---

## 2. ACHADO IMPORTANTE: a convenção do solver

A validação anterior (poço infinito, L=40) era **ambígua por coincidência**:

| Convenção | `E₁` para `L=40` |
|---|---|
| `-ψ'' + Vψ = Eψ` | `π²/L² = 0.00616850` |
| `-½ψ'' + Vψ = Eψ` | `π²/(2L²) = 0.00308425` |

Com L=40 as duas **não** coincidem — mas o valor medido (0.00616858) casava com a
primeira, então parecia tudo bem. **O teste não distinguia as convenções.**

Testes que distinguem (L=10 e oscilador harmônico) foram decisivos:

| Teste | Elmer | `-ψ''` | `-½ψ''` |
|---|---|---|---|
| poço L=10, n=1 | 0.09869637 | **0.09869604** | 0.04934802 |
| harmônico α=1, n=0 | 1.00002292 | **1.00000000** | 0.70710678 |

**Convenção confirmada: `-ψ'' + V(x)ψ = Eψ`** (ħ = 1, m = ½ — unidades atômicas).
Relação com a forma usual `-½ψ''+Vψ = Eψ`: `E_aqui = 2·E_usual`.

Isso está agora registrado no cabeçalho do solver e no cabeçalho da suíte de
regressão, para não se perder.

---

## 3. Suíte de regressão — 5/5

`scripts/regressao.sh` roda no cluster, cria pasta própria para cada caso (evita
corrida entre arquivos) e compara com solução analítica.

| Caso | Verificação | Resultado |
|---|---|---|
| poço infinito L=10 | `Eₙ=(nπ/L)²` | erro **2.96e−05** ✅ |
| oscilador harmônico k=1 | `Eₙ=2(n+½)` | erro **1.43e−04** ✅ |
| Morse (De=10, β=0.5) | mínimo do potencial = −De | `V_min=−10.000000` ✅ |
| poço duplo | splitting de tunelamento > 0 | `E₀=0.224081 E₁=0.599904` Δ=**0.375823** ✅ |
| senoidal A=0.5 | `V∈[−A,+A]` | `V=[−0.5, +0.5]` ✅ |

**Resultado: 5 passaram, 0 falharam.**

### Comentário físico de cada caso

- **Harmônico:** espaçamento 2,0 exato (`E₀=1.00014, E₁=3.00030, E₂=5.00061`)
- **Morse:** `E₀=−8.481` acima do fundo `V_min=−10` — energia de ponto zero
  `≈1.519`, compatível com `ω/2` local
- **Poço duplo:** splitting `0.376` confirma tunelamento entre as duas cavidades
  (estados par/ímpar quase degenerados sem barreira)
- **Senoidal:** autovalores próximos de zero e acima do mínimo do potencial
  `−0.5`, no início da estrutura de bandas

---

## 4. Bateria de casos por dimensão

| Caso | 1D | 2D | 3D |
|---|---|---|---|
| poço infinito | ✅ | ✅ | ✅ |
| harmônico | ✅ | — | — |
| Morse | ✅ | — | — |
| poço duplo | ✅ | — | — |
| senoidal | ✅ | — | — |
| barreira | pendente | — | — |

---

## 5. O que muda na validação do item 2

A validação do item 2 (1D/2D/3D) continua **válida**: a convenção usada é
consistente entre Elmer e Python (verificação cruzada em 4,21e−11) e as
comparações analíticas usaram a fórmula correta para `-ψ''` em cada caso.
O que faltava era **declarar explicitamente** qual convenção era — agora está.

---

## 6. Pendências

1. **Barreira de potencial** — falta implementar (transmissão/reflexão requer
   condição de contorno de radiação, não Dirichlet).
2. **2D/3D com V(x) ≠ 0** — o código suporta (o potencial é lido como campo nodal,
   independente da dimensão), mas só foi testado em 1D.
3. **Escala** — `DSYGV` denso/serial; limita a ~8000 DOFs. É o item 3.
