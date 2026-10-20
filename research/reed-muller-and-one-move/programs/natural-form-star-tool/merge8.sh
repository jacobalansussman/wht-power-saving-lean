#!/bin/bash
# merge8.sh: assemble the W = 8 starrec2 scan from its index slices (each slice file has one line per class, in
# enumeration order; process p of 2 holds the classes with index = p mod 2)
J=${1:-jobs}   # the folder where data/natural-form-star-tool/job-outputs.tar.gz was unpacked
{
  head -13000 $J/q8/remote-out/ssh2/out_q8_0.txt;  head -13000 $J/q8/remote-out/ssh2/out_q8_1.txt
  head -3000 $J/q8b/remote-out/ssh5/out_q8b_0.txt; head -3000 $J/q8b/remote-out/ssh5/out_q8b_1.txt
  grep "^S " $J/q8d/remote-out/ssh4/out_q8d_0.txt;  grep "^S " $J/q8d/remote-out/ssh4/out_q8d_1.txt
  head -3250 $J/q8c/remote-out/ssh/out_q8c_0.txt;  head -3250 $J/q8c/remote-out/ssh/out_q8c_1.txt
  head -1100 $J/q8e/remote-out/ssh5/out_q8e_0.txt; head -1100 $J/q8e/remote-out/ssh5/out_q8e_1.txt
  grep "^S " $J/q8f/remote-out/ssh/out_q8f_0.txt;   grep "^S " $J/q8f/remote-out/ssh/out_q8f_1.txt
  grep "^S " $J/q8g/remote-out/ssh4/out_q8g_0.txt;  grep "^S " $J/q8g/remote-out/ssh4/out_q8g_1.txt
} | grep "^S " > all8.txt
echo "lines $(wc -l < all8.txt) distinct classes $(awk '{print $2}' all8.txt | sort -u | wc -l)"
awk '{print $3, ($3=="P"||$3=="D") ? "" : $4}' all8.txt | sort | uniq -c
awk '$3=="P"{print $2}' all8.txt > surv8_all.txt; wc -l < surv8_all.txt
