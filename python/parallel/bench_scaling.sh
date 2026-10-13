#!/bin/bash
#==============================================================================
# Benchmark de escalabilidade MPI — autossolver paralelo (SLEPc shift-invert)
#
# Uso:  bench_scaling.sh <matriz_dir> <procs...>
# Ex.:  bench_scaling.sh ~/mtx/case 1 2 4 8
#==============================================================================
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load gnu12 openmpi4 petsc slepc 2>/dev/null

DIR=${1:-$HOME/mtx}; shift || true
PROCS=("$@"); [ ${#PROCS[@]} -eq 0 ] && PROCS=(1 2 4 8)
BIN=$HOME/bin/schrodinger_petsc
K=${K:-4}
REF=""

echo "==============================================================="
echo " Escalabilidade MPI — SLEPc shift-invert"
echo " matrices: $DIR   eigenvalues k=$K"
echo "==============================================================="
printf "  %-6s %-14s %-14s %-12s %-10s\n" "procs" "tempo_total(s)" "solucao(s)" "speedup" "efic."

for np in "${PROCS[@]}"; do
  LOG=/tmp/scale_${np}.log
  mpirun --oversubscribe -np "$np" "$BIN" -A "$DIR/A.petsc" -M "$DIR/M.petsc" \
        -k "$K" > "$LOG" 2>&1 || true
  TT=$(grep -oE "TOTAL +: +[0-9.]+ s" "$LOG" | grep -oE "[0-9.]+" | head -1)
  TS=$(grep -oE "solucao SLEPc +: +[0-9.]+ s" "$LOG" | grep -oE "[0-9.]+" | head -1)
  [ -z "$TT" ] && { echo "  np=$np FALHOU"; grep -m2 "PETSC ERROR" "$LOG"; continue; }
  if [ -z "$REF" ]; then
    REF=$TT
    printf "  %-6s %-14s %-14s %-12s %-10s\n" "$np" "$TT" "$TS" "1.00" "100%"
  else
    read -r SP EF <<<"$(awk -v r="$REF" -v t="$TT" -v n="$np" \
        'BEGIN{s=r/t; printf "%.2f %.0f%%", s, 100*s/n}')"
    printf "  %-6s %-14s %-14s %-12s %-10s\n" "$np" "$TT" "$TS" "$SP" "$EF"
  fi
done
echo "==============================================================="
