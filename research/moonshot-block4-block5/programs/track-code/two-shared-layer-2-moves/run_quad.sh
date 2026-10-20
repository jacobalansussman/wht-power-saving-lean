export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 quad.py 7 104 > out_quad_104.txt 2>&1 &
python3 quad.py 4 106 > out_quad_106.txt 2>&1 &
wait
