# Manuscript

`artigo_JSC.tex` is the manuscript submitted to the *Journal of Scientific
Computing*, with its bibliography and the compiled PDF. The LaTeX class
`sn-jnl.cls` and the style file `bst/` are the Springer Nature template,
included so that the source compiles as submitted.

## Compiling

```bash
./compila.sh
```

which runs `pdflatex`, `bibtex` and `pdflatex` twice, in that order. A full TeX
distribution is required; the template uses `cuted`, `stfloats`, `manyfoot`,
`algorithm`, `algpseudocode`, `appendix`, `threeparttable` and `wrapfig`, none
of which are in a minimal installation.

## Figures

`gera_figuras.py` regenerates the three figures in `figuras/` from the numbers
reported in the manuscript:

```bash
PYTHONPATH=/path/to/numpy-matplotlib python3 gera_figuras.py
```

## Reproducing the tables

The tables in the manuscript are reproduced by the scripts in `../verification/`
together with the case outputs in `../cases/` and `../python/casos/`.
