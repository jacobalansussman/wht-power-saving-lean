# explicit_sim.py (bound-check/rl-vprove-compute, 14 Oct 2026). Skeptic's independent check of the dilated lockstep
# tree (rl-prove PROOF.md section 3, Theorem 1G; with all ranks 1 it is ex-laser's Theorem 1 / Lemma 2).
# Written from the PROOF.md text only (no code reused from rl-prove or ex-laser).
#
# Every unit is an explicit object; every job is an explicit object with a size. Matching is done BY SIZE ONLY (the
# physical requirement), in a RANDOM order, within the class pairs the proof prescribes (round j of class s feeds class
# s+(j)). Unmatched demand is a source and must have size in (L/m, L]; unmatched supply is an error. Then a generic
# event graph (start/issue/finish of every unit) is topologically sorted (Kahn) to test deadlock-freedom independently
# of Lemma F's procedure. Phi(sources) vs Phi(exits) is compared at tau_c in 50-digit Decimal arithmetic.
# Mode "careless": the same tree, but each child's inputs are drawn by size from the pooled outputs of ALL rounds of
# its parent's whole depth (crossing class boundaries); the deadlock test must then usually find a cycle (teeth check).
import sys, random
from decimal import Decimal as D, getcontext
from math import factorial, prod
from collections import defaultdict, deque
getcontext().prec = 50

def tau_root(Wt, rounds):
    def f(t):
        return sum(D(r).ln().__mul__(t).exp() for r in rounds) - sum(w * (D(k).ln() * t).exp() for k, w in Wt.items())
    lo, hi = D(0), D(1)
    for _ in range(170):
        mid = (lo + hi) / 2
        if f(mid) > 0: lo = mid      # paid side bigger at small tau
        else: hi = mid
    return (lo + hi) / 2

def build(Wt, rounds, L, mode="lockstep", seed=1):
    rng = random.Random(seed)
    K = sorted(Wt); m = max(K); P = len(rounds)
    def c(a): return prod(k ** x for k, x in zip(K, a))
    def u(a): return factorial(sum(a)) // prod(factorial(x) for x in a) * prod(Wt[k] ** x for k, x in zip(K, a))
    def pts(rho, size):  # a in R(rho) with |a| = size
        out = []
        def rec(i, a, cc, left):
            if i == len(K) - 1:
                a2 = a + [left]; c2 = cc * K[i] ** left
                if m * rho * c2 <= L: out.append(tuple(a2))
                return
            for e in range(left + 1):
                c2 = cc * K[i] ** e
                if m * rho * c2 > L: break
                rec(i + 1, a + [e], c2, left - e)
        rec(0, [], 1, size); return out
    G = 0
    while pts(1, G + 1): G += 1
    # units: id -> dict(cls, col, slots=[(k,size)], inputs=[job], jobs=[job per round])
    units = []; jobs = []   # job: dict(size, producer=(uid,j) or None, consumer=uid or None, exit=bool)
    classes = {}            # path tuple -> list of unit ids
    def mk_class(path, rho):
        d = len(path); ids = []
        for a in pts(rho, G - d):
            for _ in range(u(a)):
                col = rho * c(a)
                units.append(dict(cls=path, col=col, slots=[(k, k * col) for k in K for _ in range(Wt[k])], inputs=[None] * sum(Wt.values()), jobs=[None] * P))
                ids.append(len(units) - 1)
        classes[path] = (ids, rho)
        return ids
    def rec_tree(path, rho):
        mk_class(path, rho)
        if len(path) < G:
            for j, r in enumerate(rounds): rec_tree(path + (j,), rho * r)
    rec_tree((), 1)
    # issue jobs
    for uid, U in enumerate(units):
        d = len(U['cls'])
        for j, r in enumerate(rounds):
            jobs.append(dict(size=r * U['col'], producer=(uid, j), consumer=None, exit=(d == G)))
            U['jobs'][j] = len(jobs) - 1
    errors = []; sources = []
    def fill(child_ids, supply_jobs):
        # match by size, random order
        pool = defaultdict(list)
        for jid in supply_jobs: pool[jobs[jid]['size']].append(jid)
        for lst in pool.values(): rng.shuffle(lst)
        slots = [(uid, i) for uid in child_ids for i in range(len(units[uid]['slots']))]
        rng.shuffle(slots)
        for uid, i in slots:
            k, size = units[uid]['slots'][i]
            if pool.get(size):
                jid = pool[size].pop(); jobs[jid]['consumer'] = uid; units[uid]['inputs'][i] = jid
            else:
                if not (L / m < size <= L): errors.append(('source size out of band', size))
                jobs.append(dict(size=size, producer=None, consumer=uid, exit=False)); units[uid]['inputs'][i] = len(jobs) - 1
                sources.append(size)
        for lst in pool.values():
            for jid in lst: errors.append(('unconsumed job', jobs[jid]['size']))
    if mode == "lockstep":
        for path, (ids, rho) in classes.items():
            if len(path) == 0: fill(ids, [])
            else:
                parent_ids, _ = classes[path[:-1]]; j = path[-1]
                fill(ids, [units[p]['jobs'][j] for p in parent_ids])
    else:  # careless: pool by depth
        bydepth = defaultdict(list)
        for path, (ids, rho) in classes.items(): bydepth[len(path)].extend(ids)
        for d in range(G + 1):
            if d == 0: fill(bydepth[0], [])
            else: fill(bydepth[d], [units[p]['jobs'][j] for p in bydepth[d - 1] for j in range(P)])
    return dict(G=G, units=units, jobs=jobs, sources=sources, errors=errors, K=K, m=m, P=P, nclasses=len(classes))

def deadlock_free(T):
    units, jobs, P = T['units'], T['jobs'], T['P']
    # events: ('S',u), ('I',u,j), ('F',u). Edges per the data dependencies.
    succ = defaultdict(list); indeg = defaultdict(int); nodes = set()
    def edge(a, b): succ[a].append(b); indeg[b] += 1; nodes.add(a); nodes.add(b)
    for uid, U in enumerate(units):
        nodes.add(('S', uid)); nodes.add(('F', uid))
        for jid in U['inputs']:
            pr = jobs[jid]['producer']
            if pr is not None: edge(('I',) + pr, ('S', uid))          # input must have been issued
        edge(('S', uid), ('I', uid, 0))
        for j in range(P):
            jid = U['jobs'][j]; nxt = ('I', uid, j + 1) if j + 1 < P else ('F', uid)
            edge(('I', uid, j), nxt)                                   # program order inside the word
            if not jobs[jid]['exit']:
                edge(('F', jobs[jid]['consumer']), nxt)                # round j's job must be COMPLETE before going on
    q = deque(n for n in nodes if indeg[n] == 0); seen = 0
    while q:
        n = q.popleft(); seen += 1
        for s in succ[n]:
            indeg[s] -= 1
            if indeg[s] == 0: q.append(s)
    return seen == len(nodes), seen, len(nodes)

def report(name, Wt, rounds, L, mode="lockstep", seed=1):
    T = build(Wt, rounds, L, mode, seed)
    tau = tau_root(Wt, rounds)
    units, jobs = T['units'], T['jobs']
    exits = defaultdict(int)
    for jb in jobs:
        if jb['exit']: exits[jb['size']] += 1
    unconsumed_nonexit = sum(1 for jb in jobs if not jb['exit'] and jb['consumer'] is None)
    phis = sum((D(s).ln() * tau).exp() for s in T['sources'])
    phie = sum(n * (D(e).ln() * tau).exp() for e, n in exits.items())
    ok, seen, tot = deadlock_free(T)
    paid = len(units) * T['P']
    print(f"{name} L={L} mode={mode} seed={seed}: tau_c={float(tau):.6f} G={T['G']} classes={T['nclasses']} units={len(units)} "
          f"paid moves issued={paid} sources={len(T['sources'])} (sizes {min(T['sources'])}..{max(T['sources'])}, L/m={L/T['m']:.2f}) "
          f"exits={sum(exits.values())} {dict(sorted(exits.items()))} errors={len(T['errors'])} {T['errors'][:2]} unconsumed={unconsumed_nonexit} "
          f"Phi_src={float(phis):.12f} Phi_exit={float(phie):.12f} |diff|={abs(phis-phie):.2E} deadlock_free={ok} ({seen}/{tot} events)", flush=True)
    return T, ok

if __name__ == "__main__":
    rules = [
        ("D rank-1 control 1+4 in 10 (Thm 1)", {3: 1, 2: 4}, [1] * 10),
        ("E rank-1 7+38 in 88 (Thm 1)", {3: 7, 2: 38}, [1] * 88),
        ("A {4:1,3:2} r=[1,1,1,1,2,2]", {4: 1, 3: 2}, [1, 1, 1, 1, 2, 2]),
        ("B {7:1,4:1} r=[1,1,1,1,2,2,2]", {7: 1, 4: 1}, [1, 1, 1, 1, 2, 2, 2]),
        ("C {5:1,4:3} r=[1]*6+[2]*2+[3]*2", {5: 1, 4: 3}, [1] * 6 + [2] * 2 + [3] * 2),
        ("H {3:2} r=[2,1,2] (word order mixed)", {3: 2}, [2, 1, 2]),
        ("J {5:2} r=[4,1,4] theta=.861", {5: 2}, [4, 1, 4]),
        ("X non-separated {3:1,2:4} r=[1]*6+[2]*2", {3: 1, 2: 4}, [1, 2, 1, 1, 2, 1, 1, 1]),
    ]
    Ls = {"D": (12, 30, 54), "E": (12,), "A": (40, 150, 400), "B": (150, 600), "C": (150, 400), "H": (300, 2000), "J": (600, 15000), "X": (30, 60)}
    for name, Wt, rl in rules:
        for L in Ls[name[0]]:
            for seed in (1, 2):
                report(name, Wt, rl, L, "lockstep", seed)
    print("--- careless (cross-class) matching: deadlock expected ---", flush=True)
    for name, Wt, rl in rules[:5]:
        L = {"D": 30, "E": 12, "A": 150, "B": 150, "C": 150}[name[0]]
        report(name, Wt, rl, L, "careless", 3)
