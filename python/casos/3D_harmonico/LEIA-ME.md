# Caso: `3D_harmonico`

| | |
|---|---|
| **Dimensão** | 3D |
| **Domínio** | `[-L/2, +L/2]^3` com `L = 12.0` |
| **Malha** | `14` elementos por direção |
| **Potencial** | `V = 0.25 · (x² + y² + ...)` (harmônico isotrópico) |
| **α** | `0.25` |

## Solução analítica
```
E_n = 2.sqrt(alpha)(n + d/2)
```

## Observação
Potencial V != 0. Degenerescencia (N+1)(N+2)/2.

## Como rodar

> ⚠ **Sempre nos nós computacionais. Nunca no headnode.**

### Opção 1 — submissão PBS (recomendado)
```bash
# a partir da pasta pai (casos/)
qsub -v CASO=3D_harmonico roda_caso.pbs
```
Saída no arquivo `3D_harmonico.o<jobid>`.

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
