gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
for s in 61 62 63; do python3 sweep_sig.py sig0 5 $s 8 100 sig_out0 > out_sig0b_$s.txt 2>&1 & done
wait
