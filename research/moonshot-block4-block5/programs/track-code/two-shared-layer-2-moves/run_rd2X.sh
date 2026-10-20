export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 rd2X.py 501 5 > out_rd2X_501.txt 2>&1 &
python3 rd2X.py 504 3 > out_rd2X_504.txt 2>&1 &
wait
