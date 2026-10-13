# Sintaxe do potencial no Elmer — e o obstáculo do MATC multidimensional

**Documento de referência.** Registra a sintaxe correta para definir `V(x,y,z)`
no SIF, incluindo a que **custa tempo descobrir**.

---

## 1. O essencial

### 1D — funciona direto
```
Body Force 1
  Potential = Variable Coordinate 1
    Real MATC "1.0*tx*tx"
End
```
Com **uma** variável, `tx` é o valor dela.

### 2D / 3D — USE VÍRGULA como separador
```
Body Force 1
  Potential = Variable Coordinate 1, Coordinate 2, Coordinate 3
    Real MATC "0.25*(tx*tx + tx(1)*tx(1) + tx(2)*tx(2))"
End
```
Com **várias** variáveis, `tx` é a **primeira** e `tx(1)`, `tx(2)`, … são as
seguintes. (Manual: *"We could have more than one variable. In that case, we use
the vector notation (example: tx(1) or tx(2))"*.)

**Verificado:** 3D harmônico, caixa L=8, malha 14³
→ `V_min = 0`, `V_max = 12.0` = `0.25·(16+16+16)` ✅
→ `E₀ = 1.558551`, degenerescência **tripla** em `E₁ = 2.591398` ✅

---

## 2. As sintaxes que NÃO funcionam (todas testadas)

| Sintaxe | Resultado |
|---|---|
| `Variable Coordinate 1 2` | `Can't find independent variable:[coordinate 1 2]` |
| `Variable Coordinate 1` + `tx(1)` | `MATC ERROR: Index out of bounds` |
| `Variable Coordinate 1 Coordinate 2` (espaço) | `Can't find independent variable` |
| duas keywords `Xc`/`Yc` com nomes distintos | `MATC ERROR: Undeclared identifier` |
| `Real MATC` multi-linha (3 linhas) | aceito mas **sem efeito** |

### Por que o espaço falha — a causa-raiz

No fonte do Elmer (`Lists.F90`, `ListParseStrToVars`), o separador é
**VÍRGULA**, não espaço:

```fortran
l1 = INDEX( str(l0:slen), ',' )     ! <-- separador e VIRGULA
```

Sem a vírgula, `Coordinate 1 Coordinate 2` é tratado como **um único nome** de
variável — que não existe.

> **Lição registrada:** o manual mostra `Velocity 1 = Variable Coordinate 2`
> (nome completo) mas **não deixa claro o separador**. A resposta estava no
> código-fonte, não na documentação.

---

## 3. Alternativa que também funciona (implementada no solver)

O solver também aceita construção **por componentes**, sem MATC multidimensional:

```
Body Force 1
  Potential X = Variable Coordinate 1
    Real MATC "0.25*tx*tx"
  Potential Y = Variable Coordinate 2
    Real MATC "0.25*tx*tx"
  Potential Z = Variable Coordinate 3
    Real MATC "0.25*tx*tx"
End
```
O solver soma `V = Vx + Vy + Vz` internamente. Útil porque cada termo é 1D e
pode vir de fontes diferentes.

**Ordem de precedência no solver:**
1. `Potential` (campo nodal — aceita as duas sintaxes acima)
2. `Potential X` + `Potential Y` + `Potential Z`
3. `V = 0`

---

## 4. Onde isso está no fonte

| Arquivo | O que tem |
|---|---|
| `Lists.F90` → `ListParseStrToVars` | parser do `Variable`, separador `,` |
| `Lists.F90` → `ListGetReal` | avalia a expressão por nó |
| `DefUtils.F90` → `GetRealValues` | variante por elemento (não usar fora do laço) |
| `EigenSolve.F90` → `ArpackEigenSolve` | **autossolver com bug** (`RESID` não inicializado) |

Cópias de referência do fonte Elmer 9.0:
`<ElmerRun>/...`
(referência **somente leitura** — cópia local em `00_planejamento/referencias/elmer/`)

---

## 5. Crédito

O fonte de referência está em `<instalação do Elmer>` (Elmer 9.0 com MUMPS,
instalado localmente) — **somente leitura**. A pasta `ElmerRun/trafo2026/`
já usava o padrão `.elmer_ref/` para consulta de API, que adotei aqui.
