# Logic, part 3: the standalone checker of commuting-matrix certificates, Hensel lifting, and the points at infinity

Section 16. `check_cert.py` uses only Python's `fractions` and `json`; each certificate in `data/two-certificate-theorem/certs/` lists its unknowns, equations, matrices and the substitution that relates the slice to the original system. The two skeptic programs read the certificates from `../../data/two-certificate-theorem/certs/`.

Outputs of the original runs: [data/two-certificate-theorem/](../../data/two-certificate-theorem/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `check_cert.py` | standalone checker of commuting-matrix certificates | `python3 check_cert.py ../../data/two-certificate-theorem/certs/*.json` | 4 s | under 100 MB | Python 3 only | VALID x6 |
| `skeptic_hensel_demo.py` | one-prime Hensel existence certificate on the Karatsuba slice | `python3 skeptic_hensel_demo.py` | under 1 s | under 100 MB | Python 3 only | identical |
| `skeptic_infinity_verify.py` | Bezout counts and the locus at infinity | `python3 skeptic_infinity_verify.py` | under 1 s | under 100 MB | Python 3 only | identical |

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
