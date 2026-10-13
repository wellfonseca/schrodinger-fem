# Caso: `3D_poco_cubico`

| | |
|---|---|
| **Dimensão** | 3D |
| **Domínio** | `[-L/2, +L/2]^3` com `L = 6.0` |
| **Malha** | `10` elementos por direção |
| **Potencial** | `V = 0` (poço) |
| **α** | `0.0` |

## Solução analítica
```
E = pi^2(nx^2+ny^2+nz^2)/L^2
```

## Observação
Degenerescencia tripla em E2=E3=E4.

## Como rodar

> ⚠ **Sempre nos nós computacionais. Nunca no headnode.**

### Opção 1 — submissão PBS (recomendado)
```bash
# a partir da pasta pai (casos/)
qsub -v CASO=3D_poco_cubico roda_caso.pbs
```
Saída no arquivo `3D_poco_cubico.o<jobid>`.

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
