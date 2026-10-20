# The price of the first touch, Conjecture FT, and the f = 2 screens of block-3 mixed units

The first-touch lemmas (Q, Q1, Q2, Q3, Q″) are written proofs (argued by hand, not checked; one AI skeptic read them
line by line). They rest on exhaustive tables, which these programs compute, and the evidence for Conjecture FT and
the f = 2 screens come from them too.

| File | What it shows |
|---|---|
| `blk3.py` | Block-3 mixed programs in one format, read three ways: exact formal programs, the hook module at f = 5, and the f = 2 test (the data splits into pieces of 8, 28 and 28 coordinates kept by every widened free step). |
| `vn2.py` | The large-f price of first touches by several steps (Lemma Q3) in the hook. |
| `srank_classes.py` | Exhaustive ranks and stable ranks of the single-step defect maps over the 14 engine and 1,680 widened classes: 8 for the cheapest two-row class, at least 20 for every other. |
| `fulltouch.py` to `fulltouch4.py` | Conjecture FT's evidence: sequences of touches that touch one two-row array (or two packed ones) completely, priced exactly in small letter groups (order 4 abelian; order 8 non-abelian in `fulltouch3.py` and `fulltouch4.py`). Every complete sequence costs at least 8 per array. |
| `twostep_quick.py` | Two-step first touches through zero divisors: the cheap first step (price 4) followed by the best second step from a pool. |
| `fam.py` | The f = 2 screen of 12-move programs (layered and kernel-relay shapes, every coefficient free, Levenberg-Marquardt from random starts, with plain 13-move controls). |
| `small11c.py` | The calibration: the 1+1 in 4 relay, which the hook forbids, stops at the same floor 0.25 and misses the same two characters. |
| `skeptic/` | The skeptic's independent rebuilds: `e7model.py` (Sp(6,2) as the rotation part of the Weyl group of E7) and `v105.py` (the defect table of Theorem M's 105-dimensional certificate over all 10,752 relevant engine steps, exact modulo two primes); `vplain.py` (the table against the plain programs); `vlift.py` (two-step touches of price exactly 1 really exist); `vhook.py` (the hook tables again); `vod.py` (the over-demand of small programs: Lemma Q″ gives 16 where Lemma Q gives 8 on a packed program). |

Commands (run from this folder, or from `skeptic/` for those files; numpy needed):

```
python3 srank_classes.py              # about 1.5 minutes
python3 twostep_quick.py              # under a minute
python3 fulltouch.py 0 2              # seed, number of batches; about a minute per batch (the letter group drawn depends on the seed)
python3 fulltouch2.py 3               # packed pairs; minutes (ran on a Colab machine); not re-run
python3 small11c.py 0 3               # seed, trials; about a minute (every trial stops at 0.25 or above)
cd skeptic && python3 v105.py         # 1 to 2 minutes
cd skeptic && python3 vplain.py ; python3 vod.py ; python3 vhook.py engine 6000 1   # seconds each (vhook: game, samples, seed)
```

`fam.py` is a library of shapes and a solver; the 86 twelve-move fits were driven by job scripts on Colab machines
and are not re-run for this revision. The saved floors and the closest program's signature are in
`../../data/first-touch-prices/` (`skeptic_out_myf2.txt`: the skeptic's re-evaluation of the nine saved programs,
the closest at 0.2500 with the missed characters 001 and 101). Theorem M's certificate itself (the 105-dimensional
representation of Sp(6,2)) is described in [../../../theory-of-theories-and-logic/](../../../theory-of-theories-and-logic/);
the author's own 105 table (`out_exh105.txt`) was computed from that certificate's saved representation, which is not
shipped here because of its size; the skeptic's model in `skeptic/` needs no file.

Comments inside the programs use the agents' track labels ("w7-mixed3" for the author, "v-w7-mixed3" for its
skeptic) and refer to their working notes; the section's Sources maps them.
