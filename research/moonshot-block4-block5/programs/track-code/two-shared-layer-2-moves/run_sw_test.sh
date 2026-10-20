gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 sweep.py priv 1 1 3 40 sw_test
