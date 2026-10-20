for s in 1 2 3; do python3 census2.py $s 400 priv > out_c2_priv_$s.txt & done
for s in 1 2 3 4; do python3 census2.py $s 400 cross > out_c2_cross_$s.txt & done
wait
