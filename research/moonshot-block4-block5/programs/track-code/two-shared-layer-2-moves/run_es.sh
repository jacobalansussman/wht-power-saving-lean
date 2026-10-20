gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 esearch_val.py > out_es_val.txt 2>&1
for s in 101 102 103 104 105 106; do python3 esearch.py 10 $s > out_es_$s.txt 2>&1 & done
wait
