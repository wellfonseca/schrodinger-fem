# Caso: `2D_poco_quadrado`

| | |
|---|---|
| **Dimensão** | 2D |
| **Domínio** | `[-L/2, +L/2]^2` com `L = 10.0` |
| **Malha** | `40` elementos por direção |
| **Potencial** | `V = 0` (poço) |
| **α** | `0.0` |

## Solução analítica
```
E = pi^2(nx^2+ny^2)/L^2
```

## Observação
Degenerescencia dupla em E2=E3 (modos 1,2 e 2,1).

## Como rodar

> ⚠ **Sempre nos nós computacionais. Nunca no headnode.**

### Opção 1 — submissão PBS (recomendado)
```bash
# a partir da pasta pai (casos/)
qsub -v CASO=2D_poco_quadrado roda_caso.pbs
```
Saída no arquivo `2D_poco_quadrado.o<jobid>`.

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
