#!/bin/sh
# Block 4, hook k = 8: check the depth-2 and depth-3 demand certificates on every one of the 16 characters alpha.
# Usage: ./run_depth.sh [PARALLEL] [WRAPPER]   (PARALLEL jobs at once, default 4; WRAPPER e.g. a memory-guard script)
# Each job is one run of `check_mj CERT ALPHA`; all 32 must print `ALPHA a: PASS`.  Then
#   f(g_2) >= 16 x 3,036 = 48,576  and  f(g_3) >= 16 x 4,270 = 68,320,   g_j = dir(e_1) ... dir(e_j).
cd "$(dirname "$0")"
PAR=${1:-4}; WRAP=${2:-}
[ -x check_mj ] || cc -O2 -o check_mj check_mj.c || exit 1
for c in g3 g2; do for a in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do echo "$c $a"; done; done |
  xargs -P "$PAR" -n 2 sh -c '$0 ./check_mj cert_m4_$1_M0.txt $2 > out_check_m4_$1_a$2.txt 2> err_check_m4_$1_a$2.txt' "${WRAP:-env}"
for c in g2 g3; do
  pass=$(cat out_check_m4_${c}_a*.txt 2>/dev/null | grep -c '^ALPHA .*: PASS')
  size=$(head -1 cert_m4_${c}_M0.txt | sed 's/.*size //')
  if [ "$pass" = 16 ]; then echo "$c: all 16 characters PASS -> f >= 16 x $size = $((16 * size))";
  else echo "$c: only $pass of 16 characters PASS -> NOT certified"; fi
done
