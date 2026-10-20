gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 openjointE.py 106 0 12 6 200 > out_jE_106_0.txt 2>&1 &
python3 openjointE.py 104 1 6 14 200 > out_jE_104_1.txt 2>&1 &
wait
