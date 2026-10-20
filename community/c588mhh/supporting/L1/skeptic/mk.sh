#!/bin/zsh
# v-host-0: compile ONE module with plain lean through locked.sh (fresh olean dir in this folder)
P=<verify>/checks/wht43/bound-check/v-host-0/proj
V=<verify>
F="$1"; M="${F%.lean}"; LOG=$P/../logs/$(echo $M | tr / _).log
PK=$V/.lake/packages
LP=$P/.lake/build/lib/lean:$V/.lake/build/lib/lean
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do LP=$LP:$PK/$p/.lake/build/lib/lean; done
mkdir -p $P/.lake/build/lib/lean/$(dirname $M)
LEAN_HEAVY=0 $V/checks/wht5/bin/locked.sh $LOG ${2:-1200} \
  /usr/bin/env LEAN_PATH=$LP LEAN_NUM_THREADS=2 /usr/bin/time -l \
  lean -DautoImplicit=false -R $P -o $P/.lake/build/lib/lean/$M.olean -i $P/.lake/build/lib/lean/$M.ilean $P/$F
