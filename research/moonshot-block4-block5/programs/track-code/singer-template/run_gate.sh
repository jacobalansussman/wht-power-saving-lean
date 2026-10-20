#!/bin/sh
# run_gate.sh W M n seed TEMPLATE -> out_gate_<T>_<W>_<M>_s<seed>.txt  (all-T relaxation gate, t_full.py)
cd "$(dirname "$0")"
exec python3 -B t_full.py "$1" "$2" "$3" "$4" "$5" > out_gate_$5_$1_$2_s$4.txt 2>&1
