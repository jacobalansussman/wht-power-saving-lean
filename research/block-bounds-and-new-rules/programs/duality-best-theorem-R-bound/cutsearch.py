# cutsearch.py -- cutting planes for the restricted LPs (pf-duality, 14 Oct 2026).  Works with any object I that has
# NF, K, m, id1(frames) -> variable id (None if the set is too big), pair_terms(d), uniform_eqs(), polymatroid_rows().
# Instances: DFZ Theorem 8.6 (char != 2) and Ingleton, letters = single frames, conditioning T = 0, 1 or 2 frames.
# Search: hill-climb from random starts at the current LP optimum (each position replaced by the best frame).
import sys, json, time, random
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from lpm import DFZ, ING, inst_row, solve

class Evaluator:
    def __init__(self, I):
        self.I = I; self.cache = {}
    def idof(self, fr):
        key = frozenset(fr)
        v = self.cache.get(key)
        if v is None:
            v = self.I.id1(key); v = -1 if v is None else v
            self.cache[key] = v
        return v
    def value(self, terms, letters, T, x):
        tot = 0.0
        for c, s in terms:
            fr = set(T)
            for l in s: fr.add(letters[l])
            v = self.idof(fr)
            if v < 0: return None
            if v: tot += c * x[v]
        return tot

def climb(E, terms, nl, x, rng, starts, tmax, passes=3):
    NF = E.I.NF; found = []
    for st in range(starts):
        letters = rng.sample(range(NF), nl)
        tsz = rng.choice([0, 1, 1, 2][:tmax + 2])
        T = tuple(rng.sample(range(NF), tsz))
        cur = E.value(terms, letters, T, x)
        if cur is None: continue
        for ps in range(passes):
            improved = False
            for pos in range(nl + 2):
                best, bestmove = cur, None
                if pos < nl:
                    for f in range(NF):
                        if f in letters: continue
                        L2 = list(letters); L2[pos] = f
                        v = E.value(terms, L2, T, x)
                        if v is not None and v < best - 1e-12: best, bestmove = v, ('L', f)
                elif pos == nl:   # change / add the first conditioning frame
                    for f in [None] + list(range(NF)):
                        T2 = tuple(t for t in T[1:]) if f is None else (f,) + tuple(T[1:])
                        if len(T2) > tmax: continue
                        v = E.value(terms, letters, T2, x)
                        if v is not None and v < best - 1e-12: best, bestmove = v, ('T', T2)
                else:
                    for f in [None] + list(range(NF)):
                        T2 = tuple(T[:1]) if f is None else tuple(T[:1]) + (f,)
                        if len(T2) > tmax or len(set(T2)) < len(T2): continue
                        v = E.value(terms, letters, T2, x)
                        if v is not None and v < best - 1e-12: best, bestmove = v, ('T', T2)
                if bestmove:
                    improved = True; cur = best
                    if bestmove[0] == 'L': letters[pos] = bestmove[1]
                    else: T = bestmove[1]
            if not improved: break
        if cur < -1e-7: found.append((cur, tuple(letters), tuple(T)))
    return found

def run(I, cuts, uniform, rounds, starts, seed, tmax, log=print, tag=""):
    rows = I.polymatroid_rows(log=log); eqs = I.uniform_eqs()
    E = Evaluator(I); rng = random.Random(seed)
    hist = []
    for rnd in range(rounds + 1):
        out = solve(I, rows, eqs, cuts, uniform=uniform, log=log)
        log(json.dumps({"tag": tag, "round": rnd, "uniform": uniform, "opt": out["opt"], "delta": out["delta"], "ncuts": len(cuts),
                        "active": out["active"][:12], "secs": round(out["secs"], 1)}))
        hist.append(out["opt"])
        if rnd == rounds: break
        added = 0; seen = set(nm for nm, _ in cuts)
        t0 = time.time()
        for kind, terms, nl in (("dfz", DFZ, 7), ("ing", ING, 4)):
            found = climb(E, terms, nl, out["x"], rng, starts, tmax)
            found.sort()
            for v, L, T in found[:80]:
                r = inst_row(I, terms, [(l,) for l in L], T)
                if not r: continue
                key = tuple(sorted(r.items()))
                nm = f"{kind}:L={list(L)},T={list(T)}"
                if key in seen: continue
                seen.add(key); cuts.append((nm, r)); added += 1
        log(json.dumps({"tag": tag, "round": rnd, "added": added, "search_secs": round(time.time() - t0, 1)}))
        if added == 0: break
    return out, cuts, hist

if __name__ == "__main__":
    from lpm import Interval, fano_letters
    m, K, uniform, rounds, starts, seed, tmax = (int(sys.argv[1]), int(sys.argv[2]), sys.argv[3] == 'uniform',
                                                 int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]))
    I = Interval(m, K)
    e = [1 << b for b in range(m)]
    cuts = []
    for nmT, T in (("none", ()), ("bottom", (0,))):
        r = inst_row(I, DFZ, fano_letters(I, 0, e[:3]), T)
        if r: cuts.append((f"dfz:points@bottom,T={nmT}", r))
    if m >= 4:
        x0 = I.IDX[frozenset({0, e[0]})]
        for nmT, T in (("none", ()), ("point", (x0,))):
            r = inst_row(I, DFZ, fano_letters(I, x0, e[1:4]), T)
            if r: cuts.append((f"dfz:lines@point,T={nmT}", r))
    out, cuts, hist = run(I, cuts, uniform, rounds, starts, seed, tmax, tag=f"m{m}K{K}")
    json.dump({"m": m, "K": K, "uniform": uniform, "history": hist, "opt": out["opt"], "delta": out["delta"],
               "active": out["active"], "cuts": [(nm, {str(k): v for k, v in r.items()}) for nm, r in cuts]},
              open(f"cuts_m{m}_K{K}_{'u' if uniform else 'a'}_{seed}.json", "w"))
