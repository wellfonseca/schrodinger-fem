#!/bin/bash
# Compila o artigo para o Journal of Scientific Computing.
#   bash compila.sh
set -u
cd "$(dirname "$0")"
# Optional: pick up a local TinyTeX/TeX Live if one is present. The script
# works with any pdflatex on the PATH.
for _d in "$HOME/.TinyTeX/bin/x86_64-linux" "$HOME/.texlive/TinyTeX/.TinyTeX/bin/x86_64-linux"; do
    [ -x "$_d/pdflatex" ] && export PATH="$_d:$PATH" && break
done
export LC_ALL=C

command -v pdflatex >/dev/null || { echo "ERRO: pdflatex nao encontrado no PATH"; exit 1; }

echo "--- 1/4 pdflatex ---"; pdflatex -interaction=nonstopmode artigo_JSC.tex > /dev/null 2>&1
echo "--- 2/4 bibtex   ---"; BSTINPUTS=.:./bst bibtex artigo_JSC > /dev/null 2>&1
echo "--- 3/4 pdflatex ---"; pdflatex -interaction=nonstopmode artigo_JSC.tex > /dev/null 2>&1
echo "--- 4/4 pdflatex ---"; pdflatex -interaction=nonstopmode artigo_JSC.tex > /tmp/jsc_last.log 2>&1

if [ ! -f artigo_JSC.pdf ]; then
  echo "ERRO: o PDF nao foi gerado. Veja /tmp/jsc_last.log"; exit 1
fi
PAG=$(pdfinfo artigo_JSC.pdf 2>/dev/null | awk '/^Pages/{print $2}')
echo
echo "OK — artigo_JSC.pdf, $PAG paginas"

UND=$(grep -icE "undefined (citation|reference)" /tmp/jsc_last.log 2>/dev/null | head -1); UND=${UND:-0}
ERR=$(grep -cE "^! " /tmp/jsc_last.log 2>/dev/null | head -1); ERR=${ERR:-0}
if [ "$UND" -gt 0 ] 2>/dev/null; then
  echo "AVISO: $UND referencia(s)/citacao(oes) nao resolvida(s)"
else
  echo "  sem referencias quebradas"
fi
if [ "$ERR" -gt 0 ] 2>/dev/null; then
  echo "AVISO: $ERR erro(s) de LaTeX:"; grep -E "^! " /tmp/jsc_last.log | head -5
else
  echo "  sem erros de LaTeX"
fi
exit 0
