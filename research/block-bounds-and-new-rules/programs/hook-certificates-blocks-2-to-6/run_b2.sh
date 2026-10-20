#!/bin/sh
# pf-logging, 14 Oct 2026: rank-2 block certificate at m = 4, f = 9 (generator about 500 MB, checker small)
cd "$(dirname "$0")"
date
/usr/bin/time -l ./gen_b2 4 8 1 2 cert_b2_m4.txt > out_gen_b2_m4.txt 2> err_gen_b2_m4.txt
date
/usr/bin/time -l ./check_b2 cert_b2_m4.txt > out_check_b2_m4.txt 2> err_check_b2_m4.txt
date
