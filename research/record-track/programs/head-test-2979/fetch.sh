#!/bin/zsh
# Read-only fetch of one file of a GitHub repository at a pinned commit (data only, never executed).
# usage: fetch.sh <owner/repo> <sha> <repo path> <download folder>
# Needs the GitHub command line tool gh (signed in). Prints "sha256  bytes  path".
set -e
repo=$1; sha=$2; p=$3; dl=$4
out=$dl/$p
mkdir -p "${out:h}"
gh api -H "Accept: application/vnd.github.raw" "repos/$repo/contents/$p?ref=$sha" > "$out"
echo "$(shasum -a 256 "$out" | cut -c1-64)  $(wc -c < "$out" | tr -d ' ')  $p"
