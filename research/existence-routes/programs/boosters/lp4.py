# e3-boost lp4.py: registers 2..5 with the certificates that do NOT rest on the Averaging Lemma at blocks >= 4
# (status proved / whp / computed / derived: A, B, C_k, distance sets, w-bounds' C, D, E, the block-3 hook and lifts,
# the hooks at blocks 4-6 (representation certificates, mod-p values)), then adding the exact block-4 lift, then all.
import certlib
from lp import rows
from lp2 import boost_min
grid = [i / 200 for i in range(1, 200)] + [0.9983]
sets = {"no frame families above block 3": [n for (n, M, P, s) in certlib.PROFILES if s in certlib.RIGOROUS and not (n.startswith("lift") and M >= 4)],
        "+ exact block-4 lift": [n for (n, M, P, s) in certlib.PROFILES if s in certlib.RIGOROUS],
        "+ projections (argued)": [n for (n, M, P, s) in certlib.PROFILES if s in certlib.RIGOROUS + ("argued",)]}
keep = list(certlib.PROFILES)
for label, names in sets.items():
    certlib.PROFILES[:] = [p for p in keep if p[0] in set(names)]
    for m in (3, 4, 5):
        R, nm = rows(m, certlib.RIGOROUS + ("argued",))
        neg = [t for t in grid if boost_min(m, t, R)[0] < -1e-9]
        print(f"[{label}] register {m}: {len(nm)} rows; boosters negative at {len(neg)}/{len(grid)}"
              + (f", from {min(neg)} to {max(neg)}" if neg else ""), flush=True)
certlib.PROFILES[:] = keep
