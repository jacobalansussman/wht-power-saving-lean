# 3 in 14 at block 5, next to plain: the zero-block test and Gröbner bases

What these programs show: RM9. `scan3.py` finds the 44 candidate sets (plain 3 in 15 minus one move, passing span and
deletion); `f2q.py` and `ztest.py` hold the field and the zero-block test; `run_z.py` runs it on all 44 in two
coordinate systems with must-fail and must-pass controls and writes `certs.json`; `mincert.py` finds the smallest
certificates; `make_gb.py` writes the local systems (Singular over Q and over a prime, msolve over two primes) for the
certificate's moves; `gb.py`, `gb2.py`, `run_gb2.py` run them for every time order and check each cofactor
certificate exactly; `t_sys.py` checks the equation generator against a closed form; `probe4.py` is the probe of the
four-array sets. The skeptic's programs are in `independent-check/` (`vscan.py`, `vcert.py`, `vgen.py`, `vmodel.py`,
`vgbrun.py` and others). The systems are in `../../data/three-in-14-groebner/systems.tar.gz` (unpacked into `gbin/`
by the work copy); msolve inputs there write powers as x^k.

Command re-run for this revision:

```
python3 run_z.py > out_run_z.txt    # 80 s: all 44 sets (both coordinate systems) excluded, controls as expected;
                                    # certs.json byte-identical to the recorded one
```

Not re-run: the Gröbner bases (they need Singular and msolve; `run_gb2.py` on a remote machine, about 0.1 s per order
on 5 moves where the system is refused), and the four-array probe.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/three-in-14-groebner/`, where the data files of `../../data/three-in-14-groebner/` sit next to the programs. Python programs need Python 3.9 or later and numpy; the Gröbner step calls Singular and msolve.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `extract_rm2.py` | DATA EXTRACTION ONLY. Uses the solver track's own field/frames (its coordinates) to turn each PASS line of the solver track/out_scan_frames3.txt into (targets as 0/1 matrices, the 14 moves as ... |
| `f2q.py` | my own small toolkit, written from the definitions (no code shared with the solver track or the skeptic). - GF(2^m) by log/exp tables, alpha a root of the same primitive polynomials (x^3+x+1, ... |
| `gb.py` | the exact algebra, theory-free and order-free. For an instance and a set of forms S (for an excluded instance: the forms of its zero-block certificate), write the polynomial system "the program's ... |
| `gb2.py` | the exact algebra, theory-free, one time order at a time (the program's own equations, as in the solver track's inst.py / tri_direct.py, restricted to a set U of moves). For an instance, a move ... |
| `independent-check/vcert.py` | [argument local: certificates whose tight pairs lie inside D_F, for the Groebner step]: for each of MY 44 sets, the zero-block certificate with the SMALLEST /D_F/ (my search over all 2^15 forms), ... |
| `independent-check/vcore.py` | my own code, written from the definitions in the general-m track 8.5-8.7; no code shared with the 3-in-14 track, the solver track or the skeptic. Field GF(32) = F2[x]/(x^5+x^2+1), alpha = x (the ... |
| `independent-check/vgbrun.py` | runs on a Colab. For MY 44 sets (vcerts_local.json: U = a 5-move zero-block certificate support whose tight pairs lie inside U), EVERY one of the 120 permutations of U (no symmetry reduction), my ... |
| `independent-check/vgen.py` | generic-m version of my tests (deletion, zero-block with local certificate, supports inside U), run on the track's must-fail controls (block-4 3 in 11, block-3 3 in 8) and their plain parents, ... |
| `independent-check/vmatch.py` | (1) convert the track's instances (its certs_min.json: its own 44 in basis [3,5,12,17,26], the solver track's 44 in the solver track's coordinates, the controls) to field elements, by finding for ... |
| `independent-check/vmodel.py` | two checks of MY equation generator (vpoly.equations), exact rational arithmetic. (1) Against the program itself: registers = W arrays + 3 scratch; a move s at form Q is the reflection I - 2 u_s ... |
| `independent-check/vplainprep.py` | must-pass control data. For each of my 44 sets: the plain 3 in 15 (deleted move restored), U+ = every copy of the certificate's directions and of the deleted direction, and every form support ... |
| `independent-check/vpoly.py` | my own sparse polynomials over Q and my own equation generator. Model : at form Q with move support D (in time order), the program needs X_D = A_D (I + G_D)^-1 B_D^T |
| `independent-check/vscan.py` | my own scan of plain 3 in 15 minus one move at block 5, in a basis-free model (field elements as directions, the polynomial basis for Sym). Symmetry: e -> 2^f e + j (mod 31) on the logs of ... |
| `independent-check/vscan_alt.py` | count plain multisets and passing sets per class representative WITHOUT the stabiliser reduction (distinct sorted multisets per target-class representative), to compare with the track's 7,076. |
| `independent-check/vscan_frames.py` | frames of every copy target (5-sets of directions whose forms add up to the target's form), my code. |
| `independent-check/vshape.py` | check the shape the track states for the 5-move certificates, on MY local certificates: the tight supports inside U are exactly the three pairs of Y, each with 2 moves and 2 odd arrays; F acts on ... |
| `independent-check/vsummary.py` | tally every Colab result file of this folder -> out_vsummary.txt. |
| `make_gb.py` | write the Groebner inputs (gb.py) for every instance in certs.json into gbin/. - Excluded instances: the subsystem of the forms in its zero-block certificate (tight witnesses + F). - Passing ... |
| `mincert.py` | for every excluded instance in certs.json, the zero-block certificate whose move set U (D_F plus the tight supports used) needs the fewest time orders (multiset permutations of U's directions, up ... |
| `p4_where.py` | WHICH TSING [maxU]: which targeted candidate of probe4.py is the file TSING (the t.sing a running probe last wrote) for? Rebuilds probe4.py's targeted candidate list for WHICH exactly as probe4.py ... |
| `probe4.py` | WHICH [maxU]: BY-PRODUCT, outside the brief. Does the local exact algebra of gb2.py see past the zero-block test on the smallest test-blind sharing shapes next to plain (four arrays)? WHICH = b4 ... |
| `run_gb2.py` | SHARD NSHARDS [--count]: the exact-algebra decision (gb2.py) for every instance in certs.json. Per instance, the move set U and the orders: - excluded instances (zero-block certificate): U = the ... |
| `run_z.py` | span, deletion, clique and zero-block tests (my code) on (1) my own 44 block-5 3-in-14 sets (out_scan3_b5.txt, my coordinates), (2) the solver track's 44 (data_rm2_44.json, its coordinates), (3) ... |
| `scan3.py` | m [W=3]: my own enumeration of 'plain minus one move' at block m . Target multisets {k_1..k_W} up to k -> 2^s k + c (mod n); each array on any frame of its target (W frames, a multiset of mW ... |
| `summarize.py` | collect res_gb2_*.jsonl (fetched into jobs/gb2/remote-out/*/) into the final table. |
| `summarize_p4.py` | coverage table of the four-array probe (probe4.py) from the fetched logs in jobs/p4/remote-out/. A run stopped by hand has finished at least N candidates, N = its last "tried N" line (probe4.py ... |
| `t_deg.py` | helper (no description in the file) |
| `t_gb2.py` | helper (no description in the file) |
| `t_gb3.py` | helper (no description in the file) |
| `t_plain.py` | helper (no description in the file) |
| `t_sys.py` | checks of gb2.system (the equation generator), exact arithmetic. (1) For random rational A, B, H and a random order, every polynomial equals the matching entry of A_D (I + H_D)^-1 B_D^T - diag(P) ... |
| `u3ctl.py` | the must-pass control U_3 = 2 in 7 (block 3, targets I, I, all seven directions, lex order 1..7) through the same exact algebra as the 44 (gb2.system with U = all 7 moves): Singular std over Q and ... |
| `z4b4.py` | my zero-block test on my own block-4 4-in-15 near-plain list (out_scan4_b4.txt). |
| `z4b5.py` | my span/deletion/clique/zero-block tests on the skeptic's nine block-5 4-in-19 survivors. |
| `ztest.py` | my own implementation of the clique test (C) and the general-m track's zero-block test 8.7 (Z), with certificates, for a multi-target configuration (targets T_l, moves z_i; any time order). Facts ... |
