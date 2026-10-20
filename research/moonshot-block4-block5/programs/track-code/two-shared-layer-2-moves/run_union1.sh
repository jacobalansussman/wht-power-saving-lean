export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 openunion.py 104 1 400 160 1500 > out_U_104_1.txt 2>&1 &
python3 openunion.py 106 0 400 160 1500 > out_U_106_0.txt 2>&1 &
wait
