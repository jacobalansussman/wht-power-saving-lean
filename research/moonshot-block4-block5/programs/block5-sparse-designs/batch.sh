#!/bin/bash
# batch.sh <list> <tag>: for each line "ztries seed args... DRV": rebuild the design (dumpsys, CAP small), then the exact
# layer-3 ball test (v4gb leaf-peeled k=2, both outputs, d=4 rows, all alpha != 0) with msolve AND Singular, char p.
export FORCE4=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=0 LMSTARTS=3 LMITERS=80 PROG=0 MODPMAX=4e8 LMMAX=3000 ZC_MAX=4
CAP=${CAP:-2000}; export CAP
P=2147483629
while read zt s drv rest; do
  ZTRIES=$zt DRV=$drv timeout 900 python3 dumpsys.py $s $rest > out_$2_s$s.txt 2>&1
  if [ -f sys_d3_$s.pkl ]; then
    INV=1 python3 v4gb.py ball sys_d4_$s.pkl sys_d3_$s.pkl 2 B${s}_k2.sing $P both >> out_$2_s$s.txt 2>&1
    python3 sing2ms.py B${s}_k2.sing B${s}_k2.ms >> out_$2_s$s.txt 2>&1
    timeout 600 msolve -t 1 -f B${s}_k2.ms -o B${s}_k2.res > /dev/null 2>&1
    timeout 600 Singular -q B${s}_k2.sing 2>&1 | grep -E 'VERDICT|time' > B${s}_k2.sout
    echo "s$s: d=4 $(grep 'd=4 round' out_$2_s$s.txt | tail -1 | cut -c1-90) | msolve $(head -c 40 B${s}_k2.res | tr '\n' ' ') | Singular $(tr '\n' ' ' < B${s}_k2.sout)" >> out_$2_summary.txt
  else
    echo "s$s: no d=3 system: $(grep -E 'd=4 round|STOP|result|no second|Traceback' out_$2_s$s.txt | tail -2 | tr '\n' ' ' | cut -c1-200)" >> out_$2_summary.txt
  fi
done < $1
