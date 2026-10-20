#!/bin/sh
H=
while read -r a; do $H python3 -W ignore mylp.py $a; done <<'L'
R1 3 eq 30 60 120 240
R1 3 eq 27 81 243
R1+R2 3 eq 30 60 120 240 243
R1+R2b 3 eq 30 60 120 240
R1+Rs 3 eq 30 60 120 240
R1+Rs 3 ge 30 60 120 240
R1 3 ge 81 240 243
R1+R2 3 ge 81 240 243
Rs 3 ge 30 60
L
