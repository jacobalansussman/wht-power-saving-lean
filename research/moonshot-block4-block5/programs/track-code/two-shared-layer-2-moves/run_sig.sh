gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
for s in 41 42 43 44; do python3 sweep_sig.py sig1 5 $s 8 60 sig_out > out_sig1_$s.txt 2>&1 & done
for s in 51 52 53; do python3 sweep_sig.py sig0 5 $s 8 60 sig_out > out_sig0_$s.txt 2>&1 & done
wait
