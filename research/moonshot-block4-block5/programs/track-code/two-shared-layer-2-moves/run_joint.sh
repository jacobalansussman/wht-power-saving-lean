gcc -O2 -shared -fPIC lib/mp.c -o lib/libmp_lin.so 2>/dev/null
python3 openjoint.py 104 1 2 8 8 200 > out_joint_104_1.txt 2>&1
