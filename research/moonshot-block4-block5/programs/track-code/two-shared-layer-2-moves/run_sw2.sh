gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
for s in 21 22 23; do python3 sweep2.py priv 6 $s 8 60 sw2_priv > out_sw2_priv_$s.txt 2>&1 & done
for s in 31 32 33 34; do python3 sweep2.py cross 5 $s 8 60 sw2_cross > out_sw2_cross_$s.txt 2>&1 & done
wait
