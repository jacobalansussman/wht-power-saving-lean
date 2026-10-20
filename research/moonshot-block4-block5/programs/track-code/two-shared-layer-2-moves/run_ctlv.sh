gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 ctl_verify.py && python3 verify.py ctl_plain8.pkl 3 4
