#!/bin/sh
# pf-logging: information only (Lemma C says every alpha-block has the same rank for odd f): rank of M_15 at m=4, f=9
cd "$(dirname "$0")"
/usr/bin/time -l ./gen_m0 4 8 15 cert_m4_M15.txt > out_gen_m4_alpha15.txt 2> err_gen_m4_alpha15.txt
