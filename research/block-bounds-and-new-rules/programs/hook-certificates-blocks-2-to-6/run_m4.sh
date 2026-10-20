#!/bin/sh
# pf-logging: block-4 dem certificate (m=4, f=9): generate the alpha=0 minor, then check it independently.
cd "$(dirname "$0")"
date
/usr/bin/time -l ./gen_m0 4 8 0 cert_m4_M0.txt > out_gen_m4.txt 2> err_gen_m4.txt
date
/usr/bin/time -l ./check_m0 cert_m4_M0.txt > out_check_m4.txt 2> err_check_m4.txt
date
