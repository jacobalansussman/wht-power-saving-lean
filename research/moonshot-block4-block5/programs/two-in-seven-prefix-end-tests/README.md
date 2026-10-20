# The 1,261 saved 2-in-7 prefixes against every end move

**What it shows.** The 2-in-7 searches of the order-3 round saved their five-move prefixes (`prefixes/all_prefixes.jsonl`, 1,261
lines). A 2 in 7 needs two end-move classes whose coefficient sums on K's frame are independent (a proved lemma).
- `sigma0.py` (Criterion 0) lists, without solving anything, the only ingredients that can reach K's frame and their
  coefficient sums. If they span at most one direction the prefix is excluded for every end move and every input basis.
- `kfD.py` (Criterion D) proves the other ingredients dead with small certificates.
- `poolsuff.py` and `poolsuff2.py` (Lemma P) handle the prefixes that already have a term on K's frame, inside the
  P-test's group of 24,576 letters.
- `sigma0_indep.py` is a second, independent implementation of Criterion 0.

**Commands.** Build the small C library that `bm.py` loads before the first run: `cc -O2 -shared -fPIC -o libmp_mac.so mp.c` on a Mac, or `cc -O2 -shared -fPIC -o libmp_lin.so mp.c` on Linux, in the folder that holds `mp.c`. (in `v3/`). Then (numpy):
- `python3 sigma0.py prefixes/all_prefixes.jsonl out.jsonl` (13 seconds on an Apple M4 Mac mini (16 GB), 45 MB): 802 prefixes have
  `dim0` 0 or 1 and are excluded;
- `python3 kfD.py IN.jsonl OUT.jsonl` and `python3 poolsuff.py IN.jsonl OUT.jsonl [RUN:IT,...]` for the rest (`p27.jsonl`,
  `p32.jsonl` are the inputs that were used; these runs took minutes to hours each and were done on Linux machines).

**Expected output.** `out_sigma0_all.jsonl` (one line per prefix), `out_summary_table.txt`, `out_p32_table.txt`,
`out_kfDx_p27.jsonl`, `out_mincert_p27.txt` (certificate sizes).

**Re-run for this revision:** `sigma0.py` on all 1,261 prefixes, on 19 October 2026: every `dim0` value matches the saved
file (198 with 0, 604 with 1, 459 with 2). Criterion D and Lemma P were not re-run.

**Open, with the first run to do:** two prefixes (anyb G3S s61 #11 and w3 STG2 #17) are excluded only against the pool
end moves. Run `poolsuff.py p32.jsonl out.jsonl anyb_G3S_s61:11,w3_STG2:17` to completion; a line
"pool-suffices" for a prefix decides it against every end move inside the P-test's group.
