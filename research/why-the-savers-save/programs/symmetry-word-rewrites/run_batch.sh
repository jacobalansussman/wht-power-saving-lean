#!/bin/sh
cd "$(dirname "$0")"
VP=${PYTHON:-python3}
WORD=../../helpers/k8-twin/words/W-k8-574-575-580-581.json.gz $VP -I scan1.py 1 7000 > scan1_pairword.log 2>&1 &
$VP -I slotperm.py $(python3 -c "print(','.join(str(x) for x in range(112,224)))") 7 > slotperm_targets.log 2>&1 &
$VP -I exhaust.py 10 300 > exhaust_L10.log 2>&1 &
wait
