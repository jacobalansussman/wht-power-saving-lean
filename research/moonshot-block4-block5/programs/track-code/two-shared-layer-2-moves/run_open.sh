gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 opencase.py 104 1 2 150 > out_open_104_1.txt 2>&1
