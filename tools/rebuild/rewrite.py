#!/usr/bin/python3 -I
"""rewrite.py: step 3 of the rebuild.  The word of pull request #193 as ONE EXPLICIT gate list (the design is
that of the outside repository, which has no explicit list of this kind; see tools/gx/ORIGIN.md).
Start from the aligned word (aligned.py; 15,682 word roles, an ordinary gcert/0 word).  Every cube has 12 pairs
of ports {p, q} whose sum p + q and difference p - q are both formed at the 2-dimensional frame span(p, q), by
two operations on FOUR port copies (two of them dead afterwards).  Replace the two operations by the in-place
map on TWO copies Ka = p, Kb = q:      Ka := Ka + Kb     (p + q)          Kb := -2 Kb + Ka     (p - q)
and delete the two other copies with the operations that created them.  1,980 pairs, 3,960 roles removed.
Everything else (frames, reuse pairs, phases, reads) is carried over with roles and operation indices renumbered.
The only new gate: a helper gate dest := ca dest + cb ctrl with ca = -2 (its undo divides by -2).

usage: python3 -I -B rewrite.py pr193      -> WORK/cache/pr193xword.pkl, WORK/cache/pr193xword.layer.pkl
"""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import defaultdict
from f2 import perp
import rbpaths

WORK = rbpaths.work()


def main(name):
    t0 = time.time()
    D = pickle.load(open(WORK + "cache/%s.aligned.pkl" % name, "rb"))
    g, wit, word, frames, pairs, newcut, prof = D["g"], D["wit"], D["word"], D["frames"], D["pairs"], D["newcut"], D["prof"]
    h, v, spans = g["h"], g["v"], wit["spans"]
    ops, coef, R = word["ops"], word["coef"], prof["R"]
    srcrole = set(word["sources"].values())
    holds = {s: x for x, s in word["sources"].items()}
    creation, nops = {}, defaultdict(int)
    for i, (a, b, x) in enumerate(ops):
        nops[a] += 1
        nops[b] += 1
        if x <= v:                                   # fan-out copy of a port
            assert coef[i] == (1, 1) and holds[b] == x and a not in holds
            holds[a], creation[a] = x, i
    groups = defaultdict(list)
    for i, (a, b, x) in enumerate(ops):
        if v < x <= newcut and len(spans[x]) == 2:
            assert frames[i] == spans[x]
            groups[frames[i]].append(i)
    pset = set(word["phase1"])
    chron = sorted(pset) + [i for i in range(len(ops)) if i not in pset]
    pos = {i: k for k, i in enumerate(chron)}
    ren, dead_ops, newop, own = {}, set(), {}, {}
    for F, ids in groups.items():
        assert len(ids) == 2 and all(i in pset for i in ids), (F, ids)
        isum, idif = sorted(ids, key=lambda i: coef[i] != (1, 1))
        assert coef[isum] == (1, 1) and coef[idif] == (1, -1)
        A1, B1, xs = ops[isum]
        A2, B2, xd = ops[idif]
        p, q = holds[A2], holds[B2]
        assert {holds[A1], holds[B1]} == {p, q} and p != q
        four = [A1, B1, A2, B2]
        assert len(set(four)) == 4
        ka = sorted((r for r in four if holds[r] == p), key=lambda r: r not in srcrole)[0]
        kb = sorted((r for r in four if holds[r] == q), key=lambda r: r not in srcrole)[0]
        for r in four:
            if r not in (ka, kb):
                assert r not in srcrole and r in creation
                dead_ops.add(creation[r])
        for r in (B1, B2):                           # the two controls are dead after their operation
            assert nops[r] == (2 if r not in srcrole else 6)
        assert A1 not in ren and A2 not in ren
        ren[A1], ren[A2] = ka, kb                    # holders of the sum and of the difference, AFTER their gate
        first, second = sorted(ids)
        newop[first] = [((ka, kb, xs), (1, 1)), ((kb, ka, xd), (-2, 1))]     # the two gates are ADJACENT
        newop[second] = []                           # (the sum holder moves on before the old second position)
        own[A1] = own[A2] = first
    dead_roles = {a for i in dead_ops for a in [ops[i][0]]}
    keep_roles = [s for s in range(R) if s not in dead_roles]
    rmap = {s: k for k, s in enumerate(keep_roles)}

    def rr(s, i=None):
        """role of the new word; a holder is renamed only for what happens after its own pair gate"""
        if s in ren and (i is None or pos[i] > pos[own[s]]):
            s = ren[s]
        assert s not in dead_roles, (s, i)
        return rmap[s]
    imap, ops2, coef2, frames2, p1 = {}, [], [], [], []
    for i, (a, b, x) in enumerate(ops):
        if i in dead_ops:
            continue
        imap[i] = len(ops2)
        if i in newop:
            for (a2, b2, x2), c2 in newop[i]:
                ops2.append((rmap[a2], rmap[b2], x2))
                coef2.append(c2)
                frames2.append(frames[i])
                if i in pset:
                    p1.append(len(ops2) - 1)
        else:
            ops2.append((rr(a, i), rr(b, i), x))
            coef2.append(coef[i])
            frames2.append(frames[i])
            if i in pset:
                p1.append(len(ops2) - 1)
    word2 = dict(word)
    word2.update(ops=ops2, coef=coef2, R=len(keep_roles),
                 phase1=sorted(p1),
                 sources={x: rmap[s] for x, s in word["sources"].items()},
                 rootroles=[rr(s) for s in word["rootroles"]],
                 selected=[dict(z, role=rr(z["role"])) for z in word["selected"]])
    pairs2 = [[rr(a), rr(b), None if t is None else imap[t]] for a, b, t in pairs]
    moved = [[i, list(U)] for i, U in enumerate(frames2) if U != perp(wit["ann"][ops2[i][2]], h)]
    prof2 = dict(prof, R=len(keep_roles))
    print("[%s] in-place pairs %d; roles %d -> %d (physical %d -> %d); operations %d -> %d; gates with ca = -2: %d; "
          "moved frames %d  (%.1fs)"
          % (name, len(groups), R, len(keep_roles), R - len(pairs), len(keep_roles) - len(pairs), len(ops), len(ops2),
             sum(1 for c in coef2 if c[0] == -2), len(moved), time.time() - t0), flush=True)
    pickle.dump(dict(g=g, prof=prof2, wit=wit, word=word2, profile={}), open(WORK + "cache/%sxword.pkl" % name, "wb"))
    pickle.dump(dict(layer=dict(frames=moved, pairs=pairs2, sinks=[])), open(WORK + "cache/%sxword.layer.pkl" % name, "wb"))


if __name__ == "__main__":
    for nm in sys.argv[1:]:
        main(nm)
