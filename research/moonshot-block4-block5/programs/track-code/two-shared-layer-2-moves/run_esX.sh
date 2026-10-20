gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
for s in 501 502 503 504; do python3 esearchX.py 6 $s > out_esX_$s.txt 2>&1 & done
wait
