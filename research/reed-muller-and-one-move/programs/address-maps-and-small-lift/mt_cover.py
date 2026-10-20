"""mt_cover.py (rm-theory): the slice's covering test (ex2-families Theorem 12, via Lemma 11 and Theorem 10) for the
MULTI-TARGET slice with one move per direction. For a fully even form F (no array odd) the block E[D_F, D_F] is
nonsingular; if every set of <= s+2 moves of D_F lies in the support of a form G with slack(G) = |D_G| - n_odd(G) <= s,
Lemma 11 makes it singular: excluded. Also scans every 4-subset of the 15 Singer (trace placement) targets at block 4
up to Singer x Frobenius symmetry: deletion, tight cliques, covering test s = 0..3.
"""
import sys, itertools
import numpy as np
from mt_lm import Setup


def cover_test(st, smax=3):
    vals = st.vals; nodd = st.phi.sum(axis=1); D = vals.sum(axis=1)
    slack = D - nodd
    supp = [frozenset(np.nonzero(vals[f])[0].tolist()) for f in range(vals.shape[0])]
    even = [f for f in range(vals.shape[0]) if D[f] > 0 and nodd[f] == 0]
    for s in range(smax + 1):
        smallG = [supp[g] for g in range(vals.shape[0]) if D[g] > 0 and slack[g] <= s and nodd[g] > 0]
        for f in even:
            if D[f] < s + 1: continue
            Sf = sorted(supp[f]); ok = True
            for k in range(1, s + 3):
                for X in itertools.combinations(Sf, k):
                    Xs = frozenset(X)
                    if not any(Xs <= G for G in smallG): ok = False; break
                if not ok: break
            if ok: return dict(excluded=True, s=s, even_form=int(f))
    return dict(excluded=False)


if __name__ == "__main__":
    m = 4
    if len(sys.argv) > 1 and sys.argv[1] == "scan":
        n = 15; seen = set(); res = {}
        for ks in itertools.combinations(range(n), 4):
            # canonical rep under k -> k + c (Singer) and k -> 2k (Frobenius on T_k = sigma^k sigma^kT)
            reps = []
            for c in range(n):
                for e in range(4):
                    reps.append(tuple(sorted(((k + c) * (2 ** e)) % n for k in ks)))
            rep = min(reps)
            if rep in seen: continue
            seen.add(rep)
            st = Setup(m, "k:" + ",".join(map(str, rep)), "all", "lex")
            et = st.exact_tests(); ct = cover_test(st)
            res[rep] = (et, ct)
            print(rep, "tight", et["tight"], "clique", et["clique"], "excluded(del/clique)", et["excluded"], "cover", ct, flush=True)
        surv = [r for r, (et, ct) in res.items() if not et["excluded"] and not ct["excluded"]]
        print(f"classes {len(res)}; surviving deletion + cliques + cover s<=3: {len(surv)}: {surv}")
    else:
        for spec in sys.argv[1:]:
            st = Setup(m, spec, "all", "lex"); print(spec, st.exact_tests(), cover_test(st))
