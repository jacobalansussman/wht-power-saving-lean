#!/bin/sh
# Zeta (subset-sum) kernel T = [[1,1],[0,1]] on 8 entries in OpenAI's finite-win game: no word with 11 or fewer
# calls computes T (x) T (x) T, so Yates's 12 additions are optimal there.
# Run from this folder:  sh run_all.sh      (needs a C compiler and a Python 3 with numpy and scipy >= 1.9;
# set PYTHON=/path/to/python3 if your default python3 lacks them).  About 15 minutes on one Apple M4 core pair,
# under 200 MB of memory.
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
cc -O2 -o tdfs tdfs.c
echo "== 1. output-order classes (mod permutations of the 3 address bits) and their cut-LP values"
$PY -I lp_pis.py                                   # writes pis_b3_lp11.txt: the 1,300 classes with LP <= 11
echo "== 2. exact integer certificates for every class with LP > 11 (they need at least 12 calls)"
$PY -I -W ignore lp_cert.py
echo "== 3. dual bounds per class (used for pruning inside the search)"
$PY -I -W ignore lp_duals.py pis_b3_lp11.txt pis_b3_lp11_y.txt
awk -F'# LP ' '{print $2 "\t" $0}' pis_b3_lp11_y.txt | sort -k1,1 -g -r | cut -f2- > pis_desc.txt
awk 'NR%2==1' pis_desc.txt > part_a.txt; awk 'NR%2==0' pis_desc.txt > part_b.txt
echo "== 4. sanity: every correct 12-call word with the identity order must survive the prunes"
grep "^0 1 2 3 4 5 6 7 " pis_b3_lp11_y.txt > pi_id_y2.txt
./tdfs 3 12 0 50 300000000 1 1 < pi_id_y2.txt > surv_g12_id.txt 2> log_g12_id.txt; tail -1 log_g12_id.txt
./tdfs 3 0 < yates_words.txt
echo "== 5. the exhaustive search: 11 calls, all 1,300 classes, two processes"
./tdfs 3 11 1 100000 300000000 1 1 < part_a.txt > surv_a.txt 2> log_a.txt &
./tdfs 3 11 1 100000 300000000 1 1 < part_b.txt > surv_b.txt 2> log_b.txt
wait
echo "classes complete: $(cat log_a.txt log_b.txt | grep -c ' complete')  aborted: $(cat log_a.txt log_b.txt | grep -c ABORTED || true)"
echo "survivors: $(cat surv_a.txt surv_b.txt | wc -l)"
tail -1 log_a.txt; tail -1 log_b.txt
echo "== 6. fewer calls (9 and 10), all classes"
for g in 9 10; do
  ./tdfs 3 $g 1 100000 0 1 1 < pis_desc.txt > surv_g$g.all.txt 2> log_g$g.all.txt
  echo "G=$g: classes $(grep -c complete log_g$g.all.txt) complete; $(tail -1 log_g$g.all.txt)"
done
echo ALL_DONE
