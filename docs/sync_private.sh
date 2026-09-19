#!/usr/bin/env bash
# Mirror the analysis to the PRIVATE repo, leaving the public fork alone.
#
# WHAT GOES PRIVATE: docs/ (at the root of the private repo) plus notebooks/,
# which sits outside docs/ in the public tree and so is missed by any
# docs-only split. Everything else -- the .f90 files, Makefile, examples --
# stays public and is not touched here.
#
# WHY NOT `git subtree push`. That recomputes the split from docs/ alone, so it
# would drop notebooks/ and then refuse the push as non-fast-forward. This
# builds the tree explicitly instead: read docs/ into a temporary index, graft
# notebooks/ alongside it, and commit that onto the existing private tip. The
# result is the same shape every time and cannot silently lose a path.
#
# The working tree and the real index are never touched -- everything happens
# in a scratch GIT_INDEX_FILE -- so this is safe to run mid-edit.
#
# Usage:   bash docs/sync_private.sh ["commit message"]
set -euo pipefail

SRC=cooper-postshutin          # branch the analysis is committed on
REMOTE=docs-private            # git@github.com:nataliaberrios/cooper-docs.git
RBRANCH=main
LOCAL=docs-only                # local mirror of the private branch
MSG=${1:-"Sync analysis from ${SRC}"}

cd "$(git rev-parse --show-toplevel)"

# Refuse to run if the private remote is missing or, worse, points at the fork.
# `git remote get-url` is git >= 2.7; this cluster has 1.8.3.1.
url=$(git config --get "remote.${REMOTE}.url" || true)
[ -n "$url" ] || { echo "no '$REMOTE' remote. add it first." >&2; exit 1; }
case "$url" in
  *cooper-docs*) ;;
  *) echo "REFUSING: '$REMOTE' is $url, which is not the private repo." >&2
     exit 1 ;;
esac

git fetch -q "$REMOTE" "$RBRANCH"
parent=$(git rev-parse FETCH_HEAD)
# Keep the local mirror pointed at the real remote tip, so a stale local
# branch cannot silently rewrite private history.
git update-ref "refs/heads/$LOCAL" "$parent"

idx=$(mktemp /tmp/syncidx.XXXXXX); rm -f "$idx"
export GIT_INDEX_FILE="$idx"
git read-tree "$SRC:docs"
git ls-tree -r "$SRC" notebooks/ \
  | awk '{printf "%s %s\t%s\n", $1, $3, $4}' \
  | git update-index --index-info
tree=$(git write-tree)
unset GIT_INDEX_FILE; rm -f "$idx"

if [ "$tree" = "$(git rev-parse "$parent^{tree}")" ]; then
  echo "private repo already matches $SRC:docs + notebooks/. nothing to do."
  exit 0
fi

commit=$(printf '%s\n' "$MSG" | git commit-tree "$tree" -p "$parent")
git update-ref "refs/heads/$LOCAL" "$commit"
git push -q "$REMOTE" "$LOCAL:$RBRANCH"

echo "pushed $(git rev-parse --short "$commit") to $REMOTE/$RBRANCH"
echo "  docs files:     $(git ls-tree -r --name-only "$commit" | grep -vc '^notebooks/')"
echo "  notebook files: $(git ls-tree -r --name-only "$commit" | grep -c '^notebooks/')"
