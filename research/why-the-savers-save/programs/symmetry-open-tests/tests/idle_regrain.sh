#!/bin/sh
# t-open test I-idle: what-if prices of the #1042 profile with the idle remainder re-grained (exact hub pricer copy).
# #1042 = 144 x local {1:2412,2:418,3:432,4:142,5:82,8:8} + idle {1:2249,2:1,3:400,7:3,16:47,48:4850}; m 73, W 14336, deficit 5120.
# Halving constraint of the moment engine: m^2 > 2 r^2, so r <= 51 at m = 73 (48 = 6 whole 8-dim stages is the largest stage-aligned child).
P=checks/wht43/bound-check/m-anatomy/pricer.py
B='"2":60193,"3":62608,"4":20448,"5":11808,"7":3,"8":1152,"16":47'
for X in '"1":349577,"48":4850' '"1":347328,"49":2249,"48":2601' '"1":349578,"49":4751' '"1":349577,"50":4656' '"1":349577,"51":4564,"36":1' '"1":349577,"51":3232,"42":1616,"48":2' '"1":349577,"51":4040,"33":808,"48":2'; do
  printf '%s  ->  ' "$X"; python3 -I $P --profile "{$B,$X}" --m 73 --W 14336 --deficit 5120 2>&1 | grep -E '"b_decimal"|kappa_decimal' | tr -d '\n '; echo
done
