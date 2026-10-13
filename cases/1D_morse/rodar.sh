#!/bin/bash
# Roda o caso. Ver ../LEIA-ME.md para pre-requisitos.
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load gnu12 openmpi4 openblas 2>/dev/null
export ELMER_HOME=$HOME/elmer_build/install_gnu
export ELMER_SOLVER_HOME=$ELMER_HOME/share/elmersolver
export LD_LIBRARY_PATH=$HOME/schrodinger_P1_2/build:$HOME/elmer_build/build_gnu/fem/src:$ELMER_HOME/lib/elmersolver:$LD_LIBRARY_PATH
ES=$HOME/elmer_build/build_gnu/fem/src/ElmerSolver
cd "$(dirname "$0")"
[ -d malha ] || python3 gerar_malha.py
$ES caso.sif > run.log 2>&1
echo "--- autovalores ---"
grep -E "SchrodingerSolver: *E\(" run.log || { echo "  FALHOU"; grep -iE "ERROR" run.log | head -5; }
