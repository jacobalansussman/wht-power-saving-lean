cd "$(dirname "$0")"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
python3 -I batch.py out_B4.jsonl specs_B4.json > log_B4.txt 2>&1
echo DONE-B4 >> log_B4.txt
