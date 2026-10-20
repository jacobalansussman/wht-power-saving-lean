#!/bin/zsh
# v-host-1 (skeptic of h-L2): clean rebuild of h-L2's XHost modules + Shea's Bank401/Adapt.lean (copied from Shea's tree)
# + an independent audit module VHost/Audit.lean (end-to-end composition). Light Lean, one module at a time via locked.sh.
B=<verify>/checks/wht43/bound-check/v-host-1; V=<verify>
L=$B/proj/.lake/build/lib/lean/Work; mkdir -p $L
for d in $V/.lake/build/lib/lean/Work/*; do b=$(basename $d); [[ "$b" == "Bank401" ]] && continue; [[ -e "$L/$b" ]] || ln -s "$d" "$L/$b"; done
cd $B/proj
for f in XHost/Core.lean XHost/Transport.lean XHost/HostBank.lean Work/Bank401/Adapt.lean XHost/Slot.lean VHost/Audit.lean; do
  $B/mk.sh $f 600 > /dev/null 2>&1 < /dev/null; ec=$?
  LG=$B/logs/$(echo ${f%.lean} | tr / _).log
  echo "$f ec=$ec $(grep -E '^EXIT' $LG | cut -d' ' -f1,2) peak=$(grep 'peak memory' $LG | awk '{print $1}') err=$(grep -c error $LG) sorry=$(grep -ci sorry $LG)"
done
grep -h "depends on axioms" $B/logs/*.log | sed 's/.*depends on axioms: //' | sort | uniq -c
