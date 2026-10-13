# Caso: `1D_harmonico` (Elmer)

| | |
|---|---|
| **Dimensão** | 1D |
| **Domínio** | `L = 20.0`, origem em `-10.0` |
| **Malha** | `400` elementos por direção |
| **Potencial** | `V = V = x^2` definido no `caso.sif` como `Variable Coordinate 1` + `Real MATC` |

## Solução analítica
```
E_n = 2(n + 1/2)
```

## Resultado esperado
E0 = 1.0. Potencial V != 0 em 1D.

## Como rodar

> ⚠ **Sempre nos nós computacionais. Nunca no headnode.**

### Opção 1 — submissão PBS (recomendado)
```bash
# a partir desta pasta
qsub roda_caso.pbs
```
Saída no arquivo `elmer_1D_harmonico.o<jobid>`.

### Opção 2 — interativo dentro de um nó
```bash
qsub -I -l select=1:ncpus=4 -q BIRA_Q3     # abre sessao no no
python3 gerar_malha.py                     # cria ./malha
$ES caso.sif > run.log 2>&1                # $ES definido em rodar.sh
grep "E(" run.log
```

## Conteúdo da pasta
| Arquivo | Função |
|---|---|
| `caso.sif` | o caso — malha, potencial e solver |
| `gerar_malha.py` | gera `./malha` no formato Elmer |
| `rodar.sh` | roda localmente (útil dentro de um nó) |
| `roda_caso.pbs` | script PBS para rodar num nó |
| `LEIA-ME.md` | este arquivo |
| **`TUTORIAL.md`** / **`TUTORIAL.pdf`** | **tutorial do caso em português** (estilo `ElmerTutorials-nonGUI`) |
| **`TUTORIAL.en.md`** / **`TUTORIAL_EN.pdf`** | **tutorial em inglês** |
| `command.txt` | todos os comandos de terminal, passo a passo |

## Pré-requisito
O binário `SchrodingerSolver.so` precisa estar instalado. Duas builds existem:

| Build | MUMPS | Caminho |
|---|---|---|
| `install_gnu` | não | `~/elmer_build/install_gnu` |
| `install_mumps` | **sim** | `~/elmer_build/install_mumps` |

O `roda_caso.pbs` usa `install_mumps` por padrão. Para usar a outra, edite
`ELMER_INSTALL` no script.

## Explicação física e validação

| Arquivo | O que é |
|---|---|
| `EXPLICACAO_FISICA.md` / `.pdf` | explicação física completa: o que o solver resolve, por que o print do ParaView parece só um traço, e as cinco checagens de validação |
| `VALIDACAO_harmonico.png` | figura pronta: potencial + nível + autofunção, numérico × analítico, erro e tabela |
| `valida_harmonico.py` | lê o `.vtu` e valida contra a solução analítica (forma, largura, normalização, virial) |
| `gera_grafico.py` | regenera a figura |

```bash
# validar
python3 valida_harmonico.py

# regerar a figura (precisa do matplotlib do projeto)
MPLCONFIGDIR=/tmp/mpl PYTHONPATH=<site-packages> \
  python3 gera_grafico.py
```

**Resumo do resultado:** estado fundamental do oscilador harmônico — uma
gaussiana `ψ₀ = e^(−x²/2)/π^(1/4)` com `E₁ = 1` na convenção do solver
(`Eₙ = 2n+1`, espaçamento constante 2). O `ψ` sai negativo porque o sinal de
uma autofunção é arbitrário. O ParaView mostra um traço fino porque a malha é
1D e a representação `Surface` não tem área — use `Plot Over Line`.
