# The best bound Theorem R can ever give (duality)

**What it shows.** Exact dual certificates for a relaxation of the best Theorem-R bound over all certificates: at
block 3 the best bound lies between 2.8 and 169/57 = 2.9649 (an 18-row certificate); at block 4, in the widened game,
between 3.154 and 3.94261; from block 4 on at most m - 0.0574. The rows are rank inequalities valid for subspace
arrangements over any field of characteristic other than 2 (submodularity, monotonicity, Ingleton-type and
Dougherty-Freiling-Zeger rows).

**Programs** (pure Python, exact fractions, no numpy for the verifiers):
- `verify3.py CERTIFICATE.json`: checks a block-3 certificate, for example
  `python3 -I verify3.py ../../data/duality/block3/certificate_u_1.json`.
- `verify_m.py CERTIFICATE.json`: checks a certificate at block m, for example
  `python3 -I verify_m.py ../../data/duality/certificates/certificate_m4_K5_u_1.json`.
- The generators (`cert3.py`, `cert_m.py`, `interval.py`, `gen_interval.py`, `dpg.py`, `dualdfz.py`, `decode_dfz.py`,
  `dfz_check.py`, `dpg_check.py`, `cert_readable.py`) need numpy and scipy, and `cert3.py` and `cert_m.py` also
  sympy. They import the linear-programme modules `lp3.py`, `lpm.py`, `cutsearch.py` and `lift_interval.py` of the
  same folder.

**Status.** The verifiers were written by the same track as the generators; no skeptic has read this track, and
nobody else has checked the certificates. An independent check is one of the open items of the section.

**Re-run for this revision:** no.
