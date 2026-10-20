export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 rd2S.py sig1 5 43 44 > out_rd2S_sig1_a.txt 2>&1 &
python3 rd2S.py sig1 5 41 42 > out_rd2S_sig1_b.txt 2>&1 &
wait
