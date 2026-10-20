gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 openjointN.py 104 1 40 4 200 > out_jointN_104_1.txt 2>&1 &
python3 openjointN.py 106 0 40 4 200 > out_jointN_106_0.txt 2>&1 &
wait
