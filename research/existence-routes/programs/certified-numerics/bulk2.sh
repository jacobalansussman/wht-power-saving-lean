#!/bin/bash
# bulk2.sh: second pass (Gaussian coefficients now supported) over the infeasible systems with no certificate yet
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1
python3 - <<'PY' > out/bulk2_files.txt
import json, glob
done = {}
for f in glob.glob('out/bulk_*.jsonl'):
    for l in open(f):
        r = json.loads(l); done[r['file']] = r
for f in sorted(done):
    r = done[f]
    if r['singular'] == 'infeasible' and not r['nulla']: print(f)
PY
for k in 0 1; do
  awk -v k=$k 'NR%2==k' out/bulk2_files.txt > out/bulk2_list_$k.txt
  ( xargs -n 10 timeout 2000 python3 bulk.py out/bulk2_$k.jsonl < out/bulk2_list_$k.txt > out/bulk2_$k.log 2>&1; echo DONE >> out/bulk2_$k.log ) &
done
wait
