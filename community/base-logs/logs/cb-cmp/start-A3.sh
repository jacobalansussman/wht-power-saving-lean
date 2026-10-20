#!/bin/bash
# A3: the official comparator on Shea's unchanged config in a tree where none of his 354 modules is compiled
# (Mathlib from the cache only), so the comparator compiles everything itself inside jail2.sh.
# landrun passes only PATH/HOME/LEAN_ABORT_ON_PANIC, so Lake would start 8 jobs; a lake adapter first on PATH
# (inside lake env, as the bigospeedrun verifier does with its lake shim) fixes LEAN_NUM_THREADS=1.
set -x
T=/content/work/cc/t-fresh; P=$T/lean/complex588; L=A3-fourier-fromnothing
s=$(date +%s)
cp -a /content/work/shea $T || exit 91
rm -rf $P/.lake/build
echo "copy wall=$(( $(date +%s)-s ))s; olean files of the repository: $(find $P/.lake/build -name "*.olean" 2>/dev/null | wc -l)"
mkdir -p /content/work/cc/fresh-bin /content/work/cc/scratch-$L
printf "#!/bin/sh\nexport LEAN_NUM_THREADS=1\nexec /content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lake \"\$@\"\n" > /content/work/cc/fresh-bin/lake
chmod +x /content/work/cc/fresh-bin/lake
chown -R 65534:65534 $P/.lake /content/work/cc/scratch-$L
cd $P
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
echo "=== $L START $(date -u) epoch=$(date +%s) config: $(sha256sum comparator/fourier.json)"
s=$(date +%s)
env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/content/work/cc/home LANG=C.UTF-8 JAIL_TMP=/content/work/cc/scratch-$L \
  /usr/bin/time -v -o /content/work/logs/cb-cmp/$L.time /content/work/cc/jail2.sh $P /bin/sh -c "exec lake env /bin/sh -c 'PATH=/content/work/cc/fresh-bin:\$PATH; export PATH; exec comparator comparator/fourier.json'"
rc=$?
echo "=== $L END rc=$rc wall=$(( $(date +%s)-s ))s $(date -u)"
grep -E "Maximum resident|Elapsed|Exit status" /content/work/logs/cb-cmp/$L.time
echo "olean files after: $(find .lake/build -name "*.olean" | wc -l); owned by someone other than 65534: $(find .lake/build ! -user 65534 | wc -l)"
echo "git status lines: $(git status --short | wc -l)"
exit $rc
