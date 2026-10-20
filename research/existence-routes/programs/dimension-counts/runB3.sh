cd "$(dirname "$0")"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
python3 -I mk_B3.py
python3 -I batch.py out_B3.jsonl specs_B3.json > log_B3.txt 2>&1
echo DONE-B3 >> log_B3.txt
