"""t_tightcode.py: block 5, canonical Singer targets, one move per direction. Check that 'some tight form' <=> 'the trace
code {(Tr(b gamma_l))_l : b} on the targets has a word of weight 2^(m-2)', and that every nonzero linear part b is the
linear part of some (m-2)-flat form. Same random draws as mt_scan5z.py."""
import sys, random
import numpy as np
from mt_lm import Setup
m = 5
st0 = Setup(m, "k:0", "all", "lex")
F = st0.F
# flats: forms of weight 2^(m-2); their 'linear part' parity pattern over the 31 targets T_k is phi; get b via parity
wt = st0.vals.sum(axis=1)
flats = np.where(wt == 2 ** (m - 2))[0]
# parity of target T_k at form f is Tr(b_f alpha^(2k)...) -- read b's pattern from the 31 canonical targets directly
stall = Setup(m, "k:" + ",".join(map(str, range(31))), "all", "lex")
pat = stall.phi                                   # forms x 31 : parity of T_k
flat_pats = set(tuple(pat[f]) for f in flats)
all_pats = set(tuple(pat[f]) for f in range(pat.shape[0]))
print("distinct nonzero target-parity patterns over all forms:", len(all_pats) - 1, "; realised by flats:", len(flat_pats))
for W in (8, 9, 10):
    rng = random.Random(5); agree = 0; N = 40
    for it in range(N):
        ks = [0] + sorted(rng.sample(range(1, 31), W - 1)) if W > 8 else [0] + sorted(rng.sample(range(1, 31), W - 1))
        st = Setup(m, "k:" + ",".join(map(str, ks)), "all", "lex")
        et = st.exact_tests()
        sub = pat[:, ks]
        has8 = any(int(sum(p[k] for k in ks)) == 2 ** (m - 2) for p in all_pats)
        agree += int(has8 == (et["tight"] > 0))
    print(f"W={W}: 'tight form exists' == 'trace code has a word of weight 8' in {agree} of {N} draws")
