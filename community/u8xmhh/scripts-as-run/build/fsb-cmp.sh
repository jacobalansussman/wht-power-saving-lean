#!/bin/bash
# fs-build: the OFFICIAL comparator (comparator d03acab, lean4export 076e8e5, nanoda 3a24072 in /content/work/cc/bin) on
# GXD.U8XMHH in my tree, positive then negative, through MY copy of c-base's jail2.sh stand-in (namespaces, read-only
# mounts except <project>/.lake and a private /tmp, uid 65534, no capabilities, no_new_privs, socket(AF_UNIX) blocked)
# -- NOT the official landrun/systemd sandbox.  k-build's kb-cmp.sh with my paths.  After the chown below, no root Lake
# in this tree.   usage: fsb-cmp.sh
F=/content/work/fs-build; P=$F/B/lean/complex588; LOGD=/content/work/logs/fs-build/cmp; mkdir -p $LOGD
N=U8XMHH; n=u8xmhh
cd $P
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
echo "tools: $(sha256sum /content/work/cc/bin/comparator /content/work/cc/bin/lean4export /content/work/cc/bin/nanoda_bin $F/cc/jail2.sh $F/cc/nounix $F/cc/fresh-bin/lake | awk '{print substr($1,1,16), $2}' | tr '\n' ' ')"
heavy() { ps -C lean -o rss= | awk '$1>2000000' | wc -l; }
for L in fsb-$N-pos fsb-$N-neg; do
  C=comparator/fourier-$n.json; [ $L = fsb-$N-neg ] && C=comparator/fourier-$n-neg.json
  mkdir -p $F/cc/scratch-$L
  chown -R 65534:65534 $P/.lake $F/cc/scratch-$L
  w=0; while [ "$(heavy)" -ge 2 ] || [ "$(free -m | awk '/Mem:/ {print $7}')" -lt 14000 ]; do sleep 10; w=$((w+10)); done
  echo "=== $L START $(date -u) (slot waited ${w}s) config: $(sha256sum $C)"; cat $C
  echo "olean files before: $(find .lake/build -name '*.olean' | wc -l)"
  s=$(date +%s)
  env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/content/work/cc/home LANG=C.UTF-8 JAIL_TMP=$F/cc/scratch-$L \
    /usr/bin/time -v -o $LOGD/$L.time $F/cc/jail2.sh $P /bin/sh -c "exec lake env /bin/sh -c 'PATH=$F/cc/fresh-bin:\$PATH; export PATH; exec comparator $C'" > $LOGD/$L.log 2>&1
  rc=$?
  echo "=== $L END rc=$rc wall=$(( $(date +%s)-s ))s $(date -u)"
  grep -E "Building|Built|accepts|okay|does not match|Const|rror" $LOGD/$L.log | cut -c1-300 | tail -10
  grep -E "Maximum resident|Elapsed|Exit status" $LOGD/$L.time
  echo "olean files after: $(find .lake/build -name '*.olean' | wc -l)"
done
echo CMPDONE
