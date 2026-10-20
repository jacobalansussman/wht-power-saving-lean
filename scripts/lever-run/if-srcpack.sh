#!/bin/bash
# i-final: deterministic tarball of EVERY source file of the final tree (no .lake), its sha256 list, and the list of
# the 700 closure modules of Work.FourierL3945.Main (from the bigospeedrun stage, whose Submission/ holds exactly the
# closure).  Read only on the tree; no lake.
set -euo pipefail
T=/content/work/i-final/T/lean/complex588; O=/content/work/i-final/srcpack; ST=/content/work/i-final/bos-L3945/stage
rm -rf $O; mkdir -p $O
cd $T
find . -path ./.lake -prune -o -type f -print | sed 's#^\./##' | LC_ALL=C sort > $O/FILES.txt
xargs -d '\n' sha256sum < $O/FILES.txt > $O/SHA256SUMS.tree
tar --sort=name --mtime='2026-10-19 00:00:00Z' --owner=0 --group=0 --numeric-owner -czf $O/l3945-tree-sources.tgz -T $O/FILES.txt
sha256sum $O/l3945-tree-sources.tgz | sed "s#$O/##" > $O/l3945-tree-sources.tgz.sha256
(cd $ST/Submission && find . -name '*.lean' | sed 's#^\./##' | LC_ALL=C sort) > $O/CLOSURE.txt
while read -r f; do [ -f "$T/$f" ] || echo "MISSING in tree: $f"; done < $O/CLOSURE.txt > $O/closure-missing.txt
grep -F -x -f $O/CLOSURE.txt $O/SHA256SUMS.tree.paths 2>/dev/null || true
awk '{print $2}' $O/SHA256SUMS.tree > $O/tree.paths
grep -F -f <(sed 's#^#  #' $O/CLOSURE.txt) $O/SHA256SUMS.tree > $O/SHA256SUMS.closure || true
echo "tree files: $(wc -l < $O/FILES.txt)  lean: $(grep -c '\.lean$' $O/FILES.txt)  closure modules: $(wc -l < $O/CLOSURE.txt)  closure hashed: $(wc -l < $O/SHA256SUMS.closure)  missing: $(wc -l < $O/closure-missing.txt)"
cat $O/l3945-tree-sources.tgz.sha256; ls -la $O
