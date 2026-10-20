export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 rd2Xb.py 5 33 > out_rd2Xc_33.txt 2>&1 &
python3 rd2Xb.py 5 34 > out_rd2Xc_34.txt 2>&1 &
wait
