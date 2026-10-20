export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 rd2.py 104 1 101 1 101 2 > out_rd2_a.txt 2>&1 &
python3 rd2.py 106 0 104 6 106 2 > out_rd2_b.txt 2>&1 &
wait
