cd "$(dirname "$0")"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
python3 -I fullness.py "$(cat specs_full4.json)" > out_full4.jsonl 2> log_F4.txt
echo DONE-F4 >> log_F4.txt
