#!/bin/bash
P=/content/work/kcomb-pkg/k-combine-u8xm
T=/content/work/k-combine/lean/complex588
L=/content/work/logs/k-combine
bash $P/scripts/build-direct.sh $P $T $L/d > $L/direct.out 2>&1
rc=$?
echo "direct rc=$rc" >> $L/direct.out
if [ $rc -eq 0 ]; then bash $P/scripts/checks.sh U8XM $T $L/d > $L/checks.out 2>&1; echo "checks rc=$?" >> $L/checks.out; fi
echo RUNDONE >> $L/direct.out
