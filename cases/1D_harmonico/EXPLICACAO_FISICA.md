# Explicação física — caso `1D_harmonico`

Resposta às três perguntas: **o que o resultado significa fisicamente**,
**o que o print mostra** e **como validar**.

---

## 1. O que o solver resolve

O `SchrodingerSolver` resolve a equação de Schrödinger **estacionária**,
unidimensional, na forma:

$$-\psi''(x) + V(x)\,\psi(x) = E\,\psi(x)$$

Neste caso o potencial é o do **oscilador harmônico**:

$$V(x) = x^2 \qquad (\alpha = 1)$$

Isso está no `caso.sif`:

```
Body Force 1
  Potential = Variable Coordinate 1
    Real MATC "1.0*tx*tx"
End
```

E a malha vai de `x = -10` a `x = +10` com 400 elementos (401 nós). As
extremidades têm `ψ = 0` (condição de contorno de Dirichlet), o que é
fisicamente correto: o potencial cresce indefinidamente, então a partícula está
**ligada** e `ψ → 0` muito antes de `x = ±10`.

---

## 2. A física do resultado

### 2.1 O estado fundamental é uma gaussiana

A solução de `−ψ'' + x²ψ = Eψ` com menor energia é

$$\psi_0(x) = \frac{1}{\pi^{1/4}}\,e^{-x^2/2}, \qquad E_0 = 1$$

Substituindo: `ψ₀'' = (x² − 1)e^{−x²/2}/π^{1/4}`, logo
`−ψ₀'' + x²ψ₀ = (−x² + 1 + x²)e^{−x²/2}/π^{1/4} = ψ₀`. Confere: `E = 1`.

**Interpretação:** `|ψ₀(x)|²` é a densidade de probabilidade de encontrar a
partícula em `x`. É uma gaussiana centrada na origem — a partícula fica
"espalhada" em torno de `x = 0`, com largura característica `σ = 1`.

### 2.2 O espectro tem espaçamento constante — a assinatura do harmônico

Para `V = αx²`, os autovalores são

$$E_n = (2n+1)\sqrt{\alpha} \qquad \Rightarrow \qquad E_n = 1,\ 3,\ 5,\ 7,\ 9,\ \dots$$

O que o solver encontrou:

| modo | `Eₙ` (solver) | analítico `2n−1` | erro |
|---|---|---|---|
| 1 | 1.000572892 | 1 | +5.7e−04 |
| 2 | 3.001197697 | 3 | +1.2e−03 |
| 3 | 5.002447062 | 5 | +2.5e−03 |
| 4 | 7.004320695 | 7 | +4.3e−03 |
| 5 | 9.006818301 | 9 | +6.8e−03 |

O **espaçamento é constante**:

```
E(2)-E(1) = 2.000625
E(3)-E(2) = 2.001249
E(4)-E(3) = 2.001874
E(5)-E(4) = 2.002498
```

Isso é *a* característica do oscilador harmônico. Em qualquer outro potencial
os níveis não ficam igualmente espaçados. É a validação mais direta e não
depende de comparar com uma fórmula — só da **estrutura** do espectro.

> **Erro residual**: os erros crescem com `n` (5.7e−04 → 6.8e−03) porque os
> estados excitados são mais estendidos e sentem mais o erro de discretização.
> É o comportamento esperado para elementos lineares.

### 2.3 Por que `ψ` sai negativo no print

O `psi` no ParaView vai de `−0.7512` a `≈ 0`. **Isso não é erro.**

Uma autofunção é definida **a menos de um sinal**: se `ψ` resolve a equação,
`−ψ` também resolve (a equação é linear e homogênea). O solver devolveu `−ψ₀`.

O que tem significado físico é `|ψ|²`, que é sempre positivo. O sinal de `ψ`
não é observável.

### 2.4 Cuidado com a convenção de unidades

O solver usa coeficiente **1** na derivada segunda. A forma usual dos livros é

$$-\frac{\hbar^2}{2m}\psi'' + \frac{1}{2}m\omega^2x^2\psi = E\psi$$

Com `ħ = m = 1`, o potencial `V = x²` corresponde a `ω = 1` e as energias
usuais são **metade** das que o solver reporta:

$$E_{\text{usual}} = \tfrac{1}{2}E_{\text{aqui}}$$

Ou seja: `E₀` aqui é `1`, mas na convenção usual é `0.5`. **As duas estão
certas** — só mudam a escala de energia. O que importa é ser consistente, e a
relação `Eₙ = 2n+1` vale na convenção do solver.

---

## 3. O que o print mostra — e por que parece só um traço

O print tem o `.vtu` correto e a física correta. O problema é
**de visualização**:

| O que você vê | Por quê |
|---|---|
| Um **traço fino** quase invisível | a malha é **1D** (401 nós sobre uma *linha*). A representação `Surface` desenha as faces das células — e uma célula de linha 1D não tem área |
| Escala de cor `−7.5e−01` a `3.1e−15` | `psi` é negativo, máximo `−0.7512` no centro e `≈0` nas bordas |
| O traço parece quase uniforme | a variação de `ψ` é suave; sem espessura não dá para ver o gradiente de cor |

**Resumindo:** o ParaView está mostrando os dados certos, mas uma malha 1D
renderizada como `Surface` é um fio de cabelo. Não há nada de errado com o
resultado.

---

## 4. Como visualizar direito no ParaView

### Opção A — `Plot Over Line` (a melhor para 1D)

Dá um **gráfico 2D** de `psi` versus `x`, que é o que você quer ver.

1. Com `run_t0001.vtu` selecionado, menu **Filters → Data Analysis → Plot Data Over Line**
   (em versões novas: **Filters → Plot → Plot Over Line**)
2. Em `Properties`, defina os pontos da linha:
   - `Point 1` = `(-10, 0, 0)`
   - `Point 2` = `(10, 0, 0)`
   - `Resolution` = `400`
3. **Apply**
4. Aparece a janela de gráfico com `psi` no eixo vertical e o comprimento de
   arco no horizontal → uma gaussiana

### Opção B — `Spreadsheet View` (ver os números)

**View → Spreadsheet View**, e escolha `psi` como `Attribute`. Você vê a tabela
`x, y, z, psi` e pode conferir contra `0.7512·exp(−x²/2)`.

### Opção C — tornar o traço visível

1. **Filters → Alphabetical → Warp By Scalar** (fator ~2) → **Apply**
2. Depois **Filters → Alphabetical → Tube** (raio ~0.05) → **Apply**
3. Agora a autofunção tem espessura e o mapa de cores fica legível

Ou, mais simples: em `Representation`, troque `Surface` por **`Points`** e
aumente `Point Size`, ou use **`Wireframe`**.

### Opção D — a figura pronta deste projeto

Já gerei a figura que responde tudo de uma vez:

```
VALIDACAO_harmonico.png
```

Ela tem o quadro clássico (potencial + nível + autofunção), a comparação
numérico × analítico, o erro absoluto e a tabela de validação. Para regerar:

```bash
PYTHONPATH=/caminho/para/site-packages \
  python3 gera_grafico.py
```

---

## 5. Como validar — cinco checagens independentes

Rode:

```bash
python3 valida_harmonico.py
```

Ele lê o `.vtu` e faz tudo. As checagens são:

### (a) Forma da autofunção — comparação com a gaussiana analítica

| grandeza | numérico | analítico | erro |
|---|---|---|---|
| `ψ(0)` | 0.7512086808 | 0.7511255445 | 8.3e−05 |
| erro máximo | — | — | 8.3e−05 |
| **erro de forma** (amplitude normalizada) | — | — | **3.2e−05** |

O erro de forma isola a *forma* da *escala*: 3.2e−05 quer dizer que o perfil
está correto a 4 dígitos.

### (b) Largura

| onde | numérico | analítico |
|---|---|---|
| `ψ = ψ(0)/e` | −1.4000 | −√2 = −1.41421 |
| `ψ = ψ(0)·e^(−1/2)` (= `−σ`) | −1.0000 | −1.00000 |

O segundo é praticamente exato (erro 1e−11) porque `x = −1` é exatamente um nó
da malha. O primeiro tem erro 0.014, que é **o mesmo que o espaçamento da
malha** (`h = 0.05`) — ou seja, é limitação de resolução, não de física.

### (c) Normalização

$$\int_{-10}^{10}|\psi|^2\,dx = 1.000208 \quad (\text{erro } 2.1\times10^{-4})$$

A autofunção está normalizada — o solver resolve o autoproblema generalizado
`Aψ = E M ψ`, que naturalmente produz autovetores `M`-ortonormais, ou seja
`∫|ψ|²dx = 1`.

### (d) Teorema do virial

Para um potencial harmônico, o teorema do virial exige `⟨T⟩ = ⟨V⟩ = E/2`.
Multiplicando a equação por `ψ` e integrando por partes:

$$\langle T\rangle = \int|\psi'|^2dx, \qquad
  \langle V\rangle = \int x^2|\psi|^2dx, \qquad
  E = \langle T\rangle + \langle V\rangle$$

| grandeza | valor | esperado |
|---|---|---|
| `⟨V⟩` | 0.5000520684 | `E/2 = 0.5002864461` |
| `⟨T⟩` | 0.5000000122 | idem |
| `⟨T⟩+⟨V⟩` | 1.0000520806 | `E₁ = 1.0005728923` |
| **`⟨T⟩/⟨V⟩`** | **0.999896** | **1.000000** |

Essa é a checagem mais forte, porque **não usa a solução analítica da
autofunção** — é uma relação de consistência interna entre o autovalor e a
autofunção. Concordância de 1e−4.

### (e) Convergência com a malha

Refine a malha e o erro deve cair como `O(h²)` para os autovalores (elementos
lineares em 1D). Edite `1D_harmonico.geo`, mude `n` e compare:

```bash
sed -i 's/^n  = 400;/n  = 800;/' 1D_harmonico.geo
gmsh 1D_harmonico.geo -2 -o 1D_harmonico.msh
ElmerGrid 14 2 1D_harmonico.msh -out mesh
ElmerSolver caso.sif
grep "E(" run.log
```

Para `n = 400` o erro em `E₁` é `5.7e−04`; dobrando `n` deve cair por ~4.

---

## 6. Um detalhe da malha que confunde

Ao ler o `.vtu` (ou `mesh.nodes`) você vai notar que a **numeração dos nós não
é monotônica**:

```
nó 1   x = -10
nó 2   x = +10      <- aqui o x "salta"
nó 3   x = -9.95
nó 4   x = -9.90
...
nó 401 x = +9.95
```

O `ElmerGrid` numera **primeiro os pontos do CAD** (os dois `Point()` do
`.geo`) e depois os nós gerados.

**Isso não é um problema:** a *conectividade* dos elementos está correta
(elemento 1 liga os nós 1 e 3, ou seja `−10 → −9.95`; elemento 400 liga 401 e 2,
`+9.95 → +10`). O ParaView usa a conectividade, então desenha certo.

Mas ao escrever um script que lê o `.vtu` você **precisa ordenar por `x`** antes
de integrar ou traçar — foi o que `valida_harmonico.py` e `gera_grafico.py`
fazem. Se não ordenar, a integral por trapézio ainda dá aproximadamente o valor
certo por acidente (as bordas têm `ψ ≈ 0`), mas o gráfico sai errado.

---

## 7. Resumo

| pergunta | resposta |
|---|---|
| **O que é fisicamente?** | estado fundamental do oscilador harmônico: gaussiana centrada em `x = 0`, `E₁ = 1` na convenção do solver |
| **Por que `ψ < 0`?** | sinal de autofunção é arbitrário; `\|ψ\|²` é o que tem significado |
| **Por que o print é um traço?** | malha 1D + representação `Surface` = fio de cabelo sem área |
| **Como ver direito?** | `Plot Over Line`, ou a figura `VALIDACAO_harmonico.png` |
| **Como validar?** | `python3 valida_harmonico.py` — forma (3e−5), normalização (2e−4), virial (1e−4) e espaçamento constante do espectro |

---

*Documento do caso `1D_harmonico` — projeto Schrodinger (Elmer + SchrodingerSolver)*
