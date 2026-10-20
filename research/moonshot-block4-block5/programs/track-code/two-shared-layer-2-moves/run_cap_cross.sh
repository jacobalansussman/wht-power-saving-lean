export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
(python3 sweep_cap.py cross 5 31 8 400 200 1200 cap_out > out_cap_cross_31.txt 2>&1; python3 sweep_cap.py cross 5 33 8 400 200 1200 cap_out > out_cap_cross_33.txt 2>&1) &
(python3 sweep_cap.py cross 5 32 8 400 200 1200 cap_out > out_cap_cross_32.txt 2>&1; python3 sweep_cap.py cross 5 34 8 400 200 1200 cap_out > out_cap_cross_34.txt 2>&1) &
wait
