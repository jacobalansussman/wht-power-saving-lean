gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 sweep2.py cross 1 5 4 60 sw2_t2
