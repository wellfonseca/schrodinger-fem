# Caso: `3D_harmonico` (Elmer)

| | |
|---|---|
| **Dimensão** | 3D |
| **Domínio** | `L = 12.0`, origem em `-6.0` |
| **Malha** | `14` elementos por direção |
| **Potencial** | `V = V = 0.25(x^2+y^2+z^2)` definido no `caso.sif` com `Variable Coordinate 1, Coordinate 2[, Coordinate 3]` + `Real MATC` |

## Solução analítica
```
E_n = 2.sqrt(a)(n+1.5)
```

## Resultado esperado
E0 ~= 1.63. Degenerescencia TRIPLA em E2=E3=E4.

## Como rodar

> ⚠ **Sempre nos nós computacionais. Nunca no headnode.**

### Opção 1 — submissão PBS (recomendado)
```bash
# a partir desta pasta
qsub roda_caso.pbs
```
Saída no arquivo `elmer_3D_harmonico.o<jobid>`.

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
