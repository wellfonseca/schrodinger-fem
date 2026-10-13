#!/bin/bash
# Publica o repositorio no GitHub e cria a release que o Zenodo arquiva.
#
# PRE-REQUISITO: o repositorio vazio tem de existir em
#   https://github.com/wellfonseca/schrodinger-fem
# Crie-o em https://github.com/new (sem README, sem .gitignore, sem licenca,
# para evitar conflito com o que ja esta commitado aqui).
#
# Uso:
#   ./PUBLICAR.sh 1.0.0
set -euo pipefail
cd "$(dirname "$0")"
VER="${1:-1.0.0}"
REMOTE="git@github.com:wellfonseca/schrodinger-fem.git"

echo "==> remoto"
git remote remove origin 2>/dev/null || true
git remote add origin "$REMOTE"

echo "==> enviando main"
git push -u origin main

echo "==> criando a tag anotada v$VER"
git tag -a "v$VER" -m "schrodinger-fem v$VER"
git push origin "v$VER"

echo
echo "Tag v$VER enviada."
echo
echo "PROXIMO PASSO (manual, no navegador):"
echo "  1. https://github.com/wellfonseca/schrodinger-fem/releases/new"
echo "     - Tag: v$VER"
echo "     - Titulo: schrodinger-fem v$VER"
echo "     - Descricao: use RELEASE_NOTES_v$VER.md"
echo "     - Publish release"
echo "  2. O Zenodo cria o DOI automaticamente em alguns minutos."
echo "  3. Verifique em https://zenodo.org/account/settings/github/"
echo "  4. Copie o DOI (concept DOI e o DOI da versao) para o artigo:"
echo "       article/artigo_JSC.tex, secao Code availability"
