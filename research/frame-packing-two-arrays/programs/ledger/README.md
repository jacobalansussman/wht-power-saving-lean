# ledger: the ledger for every unit, checked exactly

`ledger.py` is an exact engine for the engine game at blocks 1 to 4. Registers hold formal contents (sums of exact
group elements times Gaussian rationals). It runs any word, checks whether it is an exact unit, and computes the
full ledger M = W·m − 2N + X + 2D + Z + B − P (plus one per move on an empty array) in any coordinates. It also
tests whether a word can be made common-frame, or packing-free, by a change of coordinates (`find_coherent_gauge`,
`find_nopack_gauge`, `find_strict_gauge`).

**Warning.** `ledger.py` treats scalar matrices made of letters (for example −I = XZXZ) as plain numbers. In the
real model they act as (−1)^f. So `ledger.py` certifies a word only when the number of digit columns f is 1 more
than a multiple of 4. Check all four residues of f, or raw arrays at f = 1 to 6 or 8, or certify in the frame game
(`../frame-game-completion/lift.py`). `../heights/t_scalar.py` shows a word that `ledger.py` accepts and that fails at
every even f.

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/ledger/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `ledger.py` | the engine, the ledger and the coordinate tests (a module) | imported by the others | |
| `t_known.py` | the ledger on 34 named exact units at blocks 1 to 4 (plain schedules, detours, the gauge exchange, the order-3 controls, the packing gadget), in several coordinates | `python3 t_known.py 1` | 30 s |
| `t_strict.py` | strict coherence (exact labels up to scalars) of the same cases | `python3 t_strict.py` | 30 s |
| `t_strictledger.py` | the ledger with strict labels on the named cases and on random words | `python3 t_strictledger.py` | minutes |
| `t_pauliblind.py` | a unit that is common-frame up to Paulis but strictly not: the frame ledger sees no packing, the strict one does | `python3 t_pauliblind.py` | 1 s |
| `t_random.py` | the identity on 1,280 random words and 630 random exact units at blocks 2 to 4 | `python3 t_random.py 1` | minutes |
| `t_exhaust.py` | every word up to length L over a fixed alphabet: the identity on each, and every exact unit found | `python3 t_exhaust.py 2 2 4` (block 2, 692,402 runs); `python3 t_exhaust.py 3 2 4` (block 3, 2,164,802 runs) | 6 min; 23 min |
| `t_census.py` | every exact unit among words with exactly two moves at block 1 | `python3 t_census.py 1 2 2 2`, `1 2 2 3`, `1 3 2 2` | seconds to 13 min |
| `t_cutcheck.py` | the cut lemma checked directly on 20 units at blocks 2 and 3 | `python3 t_cutcheck.py` | minutes |
| `cutconst.py` | the cut constant c_m of the frame graph: exhaustive at block 2, local search at block 3 | `python3 cutconst.py 2`; `python3 cutconst.py 3 40` | 1 s; minutes |
| `cutlp.py` | the cut-condition linear program on the block-2 frame graph for 2 to 12 arrays (needs scipy) | `python3 cutlp.py` | minutes |
| `calc_needs.py` | what a 0.1 or 0.01 unit needs at each block, against the proved floors | `python3 calc_needs.py` | 1 s |

The 2.9 million runs quoted in the README are the two `t_exhaust.py` runs above (692,402 and 2,164,802). A third
enumeration (block 2 with one scratch register) was stopped after length 3. Its 320,870 runs, all without
violation, are in `t_exhaust_m2_R3_L3_partial.out`. `t_known_float.out` is the first run of `t_known.py` with
floating-point group arithmetic, kept to show that the exact version gives the same lines. `rerun_openai_labels_h7.out`
is OpenAI's network checked with `../packing-measurement/labels_check.py` at h = 7.
