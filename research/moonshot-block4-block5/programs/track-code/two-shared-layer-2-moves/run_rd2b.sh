export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 rd2batch.py 12 111 112 > out_rd2b_A.txt 2>&1 &
python3 rd2batch.py 12 113 114 > out_rd2b_B.txt 2>&1 &
wait
