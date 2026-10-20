# Block 5, the 2-in-8 shape: designs, the corrected forced-zero closure, and exact decisions

**What it holds.** The design and decision code of the last block-5 round.
- `sp2.py`, `mwlib.py`, `endsolve.py`, `split.py`, `jdesign.py`: the frame geometry of block 5 (Sp(10,2)), the split
  trunk and the orbit designs.
- `peel2.py`: the forced-zero closure. A row with a single live producer forces that producer to zero; the closure
  spreads these zeros and audits itself for completeness. `peel.py` is the earlier version, kept as the record: it
  missed forced zeros in two ways (a column killed before it was built stayed listed as a producer; a row met while all
  its built producers were dead was dropped for good). It never produced a false exclusion, but it kept dead designs
  looking alive. Use `peel2.py`.
- `alpha3.py`, `v4torus.py`, `modsolve.py`, `sparsemod.py`: the layer-4 solution set V4 (along each orbit the injection
  coefficients form a geometric sequence whose ratio is a root of unity) and layer-3 tests on it, over finite fields.
- `v4gb.py`, `ballm.py`, `localize.py`, `dumpsys.py`, `sing2ms.py`: the one-shot ball test (rows within r steps of the
  seed, with an inverse for every injection coefficient) written for Singular and msolve. `jd_union.py`, `jd_zcover.py`:
  the union and z-cover designs. `certify.py`: checks forced-zero chains. `batch.sh`, `u.list`, `z.list`, `t11.sh`,
  `t34.sh`, `t35.sh`, `t36.sh`, `t42_1401.sh`, `t42_1402.sh`, `t46.sh`: how the runs were launched (the `t3*` and `t42*`
  scripts write the layered systems that `../block5-ball-test/vball.py` reads).
- `skeptic/vcert.py` (brute-force check of every kill chain) and `skeptic/vlm.py` (the skeptic's own numerical solver).
- `outputs/`: summary lines of the decisive runs (`t11/` holds the re-decision of the 20 old near misses).

**Needs** numpy, python-flint (for `modsolve.py`), sympy (for `allparam.py`), and Singular and msolve for the ball tests.

**Warnings.** `sing2ms.py` corrupts any entry that is a sum of two terms (no such entry occurs in this round's systems).
msolve 0.6.5 reads a product written `x*x` as `x`; write powers as `x^k` and cross-check every msolve verdict with
Singular.

**Re-run for this revision:** no. The runs took minutes to hours on Linux machines of Google's Colab service, and
python-flint, Singular and msolve are not installed on the Mac. The decisive systems are in `../block5-ball-test/`,
regenerated and checked byte for byte.
