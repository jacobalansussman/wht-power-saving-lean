export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 sweep_cap.py sig1 5 43 8 400 200 1200 cap_out > out_cap_sig1_43.txt 2>&1 &
python3 sweep_cap.py sig1 5 44 8 400 200 1200 cap_out > out_cap_sig1_44.txt 2>&1 &
wait
