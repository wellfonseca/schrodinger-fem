# Caso: `1D_poco_duplo` (Elmer)

| | |
|---|---|
| **Dimensão** | 1D |
| **Domínio** | `L = 20.0`, origem em `-10.0` |
| **Malha** | `400` elementos por direção |
| **Potencial** | `V = V = 0.01(x^2-4)^2` definido no `caso.sif` como `Variable Coordinate 1` + `Real MATC` |

## Solução analítica
```
splitting de tunelamento
```

## Resultado esperado
E0=0.224 E1=0.600. Splitting = 0.376.

## Como rodar

> ⚠ **Sempre nos nós computacionais. Nunca no headnode.**

### Opção 1 — submissão PBS (recomendado)
```bash
# a partir desta pasta
qsub roda_caso.pbs
```
Saída no arquivo `elmer_1D_poco_duplo.o<jobid>`.

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
