#!/bin/bash
ulimit -v 12000000
python3 -B -c "import chk_bf7d; print('import ok')"
python3 -B chk_bf7d.py 0 72 20000 20000 4 400000 > out_chk_s72t0.txt 2>&1 &
python3 -B chk_bf7d.py 0 71 20000 20000 4 400000 > out_chk_s71t0.txt 2>&1 &
wait
