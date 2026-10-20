gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 sweep_sig.py sig1 1 41 6 60 sig_test
