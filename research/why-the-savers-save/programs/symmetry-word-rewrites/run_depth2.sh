#!/bin/sh
cd "$(dirname "$0")"
VP=${PYTHON:-python3}
$VP -I depth2.py 10 > depth2_k8.log 2>&1 &
WORD=../../helpers/k8-twin/words/W-k8-574-575-580-581.json.gz $VP -I depth2.py 10 > depth2_pairword.log 2>&1 &
wait
