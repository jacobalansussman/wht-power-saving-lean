"""t_fastpath.py: the evaluators' fast Pauli path (index permutation and sign) against the plain tensor-power product
with the explicit matrix X^a Z^b, for every Pauli at blocks 2 and 3 and f = 1..4, and the scalar rule [i^j g] -> i^(jf)."""
import sys, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from complete import *
for m in (2, 3):
    E = Eng(m); worst = 0
    for f in range(1, 5):
        R = Raw(E, f); rn = np.random.default_rng(f)
        v = rn.normal(size=R.N) + 1j * rn.normal(size=R.N)
        for u in range(E.NL):
            worst = max(worst, np.abs(R.pauli(u, v) - R.mat(E.P[u], v)).max())
        for j in range(4):
            worst = max(worst, np.abs(R.nfree({(5 % E.NL, j): 1.0}, v) - R.mat((1j ** j) * E.P[5 % E.NL], v)).max())
    print(f"m={m}: fast Pauli path vs explicit tensor power, all {E.NL} Paulis, f = 1..4: max difference {worst:.1e}")
