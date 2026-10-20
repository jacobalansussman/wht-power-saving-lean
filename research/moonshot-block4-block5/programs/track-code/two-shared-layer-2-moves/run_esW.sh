gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 esearchW.py 2 120 201 40 > out_esW2_201.txt 2>&1 &
python3 esearchW.py 3 160 301 25 > out_esW3_301.txt 2>&1 &
python3 esearchW.py 3 160 302 25 > out_esW3_302.txt 2>&1 &
python3 esearchW.py 4 260 401 6 > out_esW4_401.txt 2>&1 &
wait
