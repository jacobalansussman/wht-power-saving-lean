#!/bin/sh
H=
while read -r a; do $H python3 -W ignore mylp.py $a; done <<'L'
R1+Rs 3 ge 27 81 243
R1+R2b 3 ge 27 81 243
R1+Rs 1 ge 27 81 243
R1+Rs 2 ge 81 243
R1+Rs 6 ge 27 81
R1+R2b 1 ge 27 81 243
R1+R2 1 ge 81 240 243
GT 1 eq 60 120
GR 1 eq 60 120
GT 1 ge 60 120
GR 1 ge 60 120
GM 1 eq 60 120
GM 1 ge 60 120
L
