gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 opencase.py 106 0 1 150 > out_open_106_0.txt 2>&1
