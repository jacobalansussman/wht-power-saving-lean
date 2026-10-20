export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 rd2Xb.py 6 601 602 > out_rd2Xb_A.txt 2>&1 &
python3 rd2Xb.py 6 603 604 > out_rd2Xb_B.txt 2>&1 &
wait
