"""lemmaR1_check.py -- o3-shapeR: computational check of Lemma R1 (PROOF.md): if on input l the third move u3 has
exactly one depth-2 term, at least one depth-1 (raw) frame other than N_l = frame(u_l), and every frame of u3 on input l
carries exactly one term, then eps_l(M4) = eps_l(M5') = 0 for every target pair.  Builds prefixes of the families
RAWR / RAWS / MIX / SECR / COS (run_fam.py's builder) and checks every input that meets the hypotheses."""
import sys, os, collections
sys.argv = [sys.argv[0], sys.argv[1], "0", sys.argv[2], "/dev/null"] + sys.argv[3:]
fam, seed, NP = sys.argv[1], int(sys.argv[3]), int(sys.argv[5]) if len(sys.argv) > 5 else 100
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_fam.py")).read().split("stats = collections.Counter()")[0])
NP = int(sys.argv[5]) if len(sys.argv) > 5 else 100
st = collections.Counter()
for k in range(NP):
    b = build()
    if b is None or isinstance(b, str):
        continue
    gens, info = b
    u3 = gens[4]; u = [gens[2], gens[3]]
    res, n1, n2 = test_prefix(W, gens)
    # Lemma R2 hypothesis: no depth-1 frame of u3 other than N_l on either input, one term per frame
    r2 = True
    for l in (0, 1):
        Nl = fk(u[l][0][1]); fr = collections.Counter(fk(M) for (ll, M, q) in u3 if ll == l)
        for (ll, M, q) in u3:
            if ll == l and W.st(M)[0] == 1 and fk(M) != Nl:
                r2 = False
        if any(v > 1 for f, v in fr.items() if f != Nl and f != W.L0):
            r2 = False
    st[(fam, "R2hyp" if r2 else "noR2", "dim eps(M4) max", max(x["d4"][0] for x in res), "dim eps(M5') max", max(x["d5"][0] for x in res))] += 1
    for l in (0, 1):
        Nl = fk(u[l][0][1])
        frames = collections.Counter(); dep = {}
        for (ll, M, q) in u3:
            if ll != l: continue
            f = fk(M); frames[f] += 1; dep[f] = W.st(M)[0]
        d2 = [f for f in frames if dep[f] == 2]; d1 = [f for f in frames if dep[f] == 1 and f != Nl]
        single = all(v == 1 for f, v in frames.items() if f != Nl and f != W.L0)
        hyp = len(d2) == 1 and len(d1) >= 1 and single
        e4 = max(x["d4"][1 + l] for x in res); e5 = max(x["d5"][1 + l] for x in res)
        # Lemma R1'' hypothesis: some single-term depth-1 frame (!= N_l) adjacent to <= 1 depth-2 frame of u3 on input l,
        # every depth-2 frame of u3 on input l single-term
        d2single = all(frames[f] == 1 for f in d2)
        r1pp = False
        for f in d1:
            if frames[f] != 1: continue
            M = next(MM for (ll, MM, qq) in u3 if ll == l and fk(MM) == f)
            nbf = set(W.nbrs(M)[1])
            if sum(1 for f2 in d2 if f2 in nbf) <= 1:
                r1pp = True
        if r1pp and d2single:
            st[(fam, "R1pp-hyp", "eps_l(M4) max", e4)] += 1
        st[(fam, "hyp" if hyp else "nohyp", "eps_l(M4) max", e4, "eps_l(M5') max", e5)] += 1
for k, v in sorted(st.items(), key=lambda kv: str(kv[0])):
    print(k, v)
