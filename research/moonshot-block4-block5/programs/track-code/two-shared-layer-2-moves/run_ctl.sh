gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 sp_check2.py
python3 ctl.py
