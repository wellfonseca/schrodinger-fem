# Caso: `1D_poco_infinito`

| | |
|---|---|
| **Dimensão** | 1D |
| **Domínio** | `[-L/2, +L/2]^1` com `L = 40.0` |
| **Malha** | `250` elementos por direção |
| **Potencial** | `V = 0` (poço) |
| **α** | `0.0` |

## Solução analítica
```
E_n = (n.pi/L)^2
```

## Observação
Valida a formulacao basica. Erro esperado < 1e-4.

## Como rodar

> ⚠ **Sempre nos nós computacionais. Nunca no headnode.**

### Opção 1 — submissão PBS (recomendado)
```bash
# a partir da pasta pai (casos/)
qsub -v CASO=1D_poco_infinito roda_caso.pbs
```
Saída no arquivo `1D_poco_infinito.o<jobid>`.

### Opção 2 — interativo dentro de um nó
```bash
qsub -I -l select=1:ncpus=4 -q BIRA_Q3     # abre sessao no no
cd $(dirname $0)
python3 run.py
```

## Conteúdo da pasta
| Arquivo | Função |
|---|---|
| `run.py` | o caso — monta, resolve e compara com o analítico |
| `roda_caso.pbs` | script PBS para rodar este caso num nó |
| `LEIA-ME.md` | este arquivo |

O ponto de entrada `run.py` importa `_comum.py` (na pasta pai) e
`../scripts/fem_schrodinger.py` (o solver compartilhado).
