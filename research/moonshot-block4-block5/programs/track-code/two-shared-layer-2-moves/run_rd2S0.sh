export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 rd2S.py sig0 5 61 62 > out_rd2S_sig0_a.txt 2>&1 &
python3 rd2S.py sig0 5 63 > out_rd2S_sig0_b.txt 2>&1 &
wait
