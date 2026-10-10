"""compile.py: from the signed addition circuit to the helper word.

  closure(g, arcs)      : frames of all nodes (annihilators = full backward intersections along uses and
        carrier arcs), the operation order, the role count R = additions + root uses - arcs, the block ledger H.
  select(g, prof, wit)  : role allocation (the word M: gates (dest, control, node) with signs), the two
        phases, and the entrance gauges with the targets they are charged to.
ADAPTED from two scripts of the outside repository (github.com/CrocSwap/integer-mult-bounds, Apache-2.0) at
commit 4a3c769.  They were read as text and never run.
  closure()  follows compile_closure() of scripts/paired_cube/closure.py statement by statement.  Notice of that
             file: "Copyright 2026 DaysSky (PR162, research/extended-links/cxlinks.py), Apache-2.0; moved here
             unchanged apart from this header and the imports (eumemic, Claude assistance)."  The file says that
             its function is compile_graph() of scripts/paired_cube/frames.py for an arbitrary legal arc list.
  select()   follows the first part of select() of scripts/paired_cube/gauges.py statement by statement.  Notice
             of that file and of frames.py: "Copyright 2026 icekylinx. Apache-2.0." and "Developed with GPT-6
             Astra assistance; integrated with Codex assistance."
Changed here: one statement to a line; some names (basis -> rref, contained -> inside, ell -> loss); 0x7fffffff
written as ROOTBIT - 1; the option repository_order and the parts that fill the outside records are left out;
the results are returned as the dictionaries prof and wit.  NOTICE, section 6, lists all files of this kind.
The indices (1-based nodes, use code 2*node+operand, root use (1<<31)|root) are those of the outside
repository, because its frozen data files (matching arcs, moved frames, reuse pairs) are indexed that way.
"""
import heapq, math
from collections import Counter
from f2 import rref, perp, inside

ROOTBIT = 1 << 31


def shift(g):
    """1-based copy: node 0 reserved"""
    args = [None] + [None if a is None else [a[0] + 1, a[1] + 1] for a in g["args"]]
    roots = [dict(r, node=r["node"] + 1) for r in g["roots"]]
    return args, roots


def closure(g, arcs):
    h, inputs = g["h"], g["inputs"]
    args, roots = shift(g)
    v, n, q = len(inputs), len(args), len(roots)
    spans = [()] * n
    for i, u in enumerate(inputs):
        spans[i + 1] = (u,)
    for x in range(v + 1, n):
        a, b = args[x]
        assert 0 < a < x and 0 < b < x
        spans[x] = rref(spans[a] + spans[b])
    rframe, rann, Y, targetH, loss = [], [], [()] * v, Counter(), 0
    for r in roots:
        x = r["node"]
        if r["kind"] == "center":
            U = spans[x]
            A = perp(U, h)
            loss += len(U)
        else:
            A = rref(inputs[t] for t in r["targets"])
            U = perp(A, h)
            for t in r["targets"]:
                assert inside(Y[t], U), "target retreat"
                targetH[len(U) - len(Y[t])] += 1
                Y[t] = U
        assert inside(spans[x], U), "root incompatible"
        rframe.append(U)
        rann.append(A)
    for t in range(v):
        assert inside(Y[t], perp((inputs[t],), h))
        targetH[h - 1 - len(Y[t])] += 1
    active = set(range(1, v + 1))
    todo = [r["node"] for r in roots]
    while todo:
        x = todo.pop()
        if x not in active:
            active.add(x)
            if args[x]:
                todo.extend(args[x])
    arcs = {int(x): int(c) for x, c in arcs}

    def use_node(c):
        return None if c >> 31 else c // 2

    def use_value(c):
        return roots[c & (ROOTBIT - 1)]["node"] if c >> 31 else args[c // 2][c & 1]
    succ, direct, indeg = [[] for _ in args], [[] for _ in args], [0] * n
    for x in sorted(active):
        if args[x]:
            for y in args[x]:
                succ[y].append(x)
                indeg[x] += 1
    for j, r in enumerate(roots):
        direct[r["node"]].extend(rann[j])
    assert len(set(arcs.values())) == len(arcs), "one donor per use"
    for x, c in arcs.items():
        assert args[x] and use_value(c) in args[x] and c not in (2 * x, 2 * x + 1)
        t = use_node(c)
        if t is None:
            direct[x].extend(rann[c & (ROOTBIT - 1)])
        else:
            succ[x].append(t)
            indeg[t] += 1
    heap = [x for x in sorted(active) if indeg[x] == 0]
    heapq.heapify(heap)
    topo = []
    while heap:
        x = heapq.heappop(heap)
        topo.append(x)
        for y in succ[x]:
            indeg[y] -= 1
            if indeg[y] == 0:
                heapq.heappush(heap, y)
    assert len(topo) == len(active), "dependency graph cyclic"
    ann, rank = [None] * n, [0] * n
    for x in reversed(topo):
        ann[x] = rref(direct[x] + [z for y in succ[x] for z in ann[y]])
        rank[x] = h - len(ann[x])
        assert inside(spans[x], perp(ann[x], h)), "value span outside frame"
    index = {x: i for i, x in enumerate(topo)}
    order = sorted(active, key=lambda x: (rank[x], index[x]))     # frame rank, then generalized topological index
    pos = {x: i for i, x in enumerate(order)}
    for x in order:
        for y in succ[x]:
            assert pos[y] > pos[x], "order is not a linear extension"
    uses = [[] for _ in args]
    for x in order:
        if args[x]:
            for j, y in enumerate(args[x]):
                uses[y].append(2 * x + j)
    for j, r in enumerate(roots):
        uses[r["node"]].append(ROOTBIT | j)
    donors = [x for x in order if args[x]]
    H = Counter()
    for x in order:
        r, deg = rank[x], len(uses[x])
        assert deg > 0
        H[r] += deg - 1
        if args[x]:
            H[h - r] += 1
            for y in args[x]:
                assert inside(ann[x], ann[y])
                H[r - rank[y]] += 1
        else:
            H[1] += 1
            H[r - 1] += 1
    for j, rt in enumerate(roots):
        r, top = rank[rt["node"]], len(rframe[j])
        assert r <= top
        if rt["kind"] == "center":
            assert r == top
            H[r] += 1
            H[h - r] += 1
        else:
            H[top - r] += 1
            H[h - top] += 1
    for x, c in arcs.items():
        val, t = use_value(c), use_node(c)
        rv, rd = rank[val], rank[x]
        top = rank[t] if t is not None else len(rframe[c & (ROOTBIT - 1)])
        assert top >= rd >= rv
        H[h - rd] -= 1
        H[rv] -= 1
        H[top - rv] -= 1
        H[top - rd] += 1
    assert min(H.values()) >= 0
    R = len(donors) + q - len(arcs)
    assert sum(r * c for r, c in H.items()) == h * R + loss
    prof = dict(h=h, v=v, R=R, q=q, c=len(donors), matched=len(arcs), loss=loss, H=H, targetH=targetH,
                source=Counter({2: v, h - 4: v, 1: v}))
    wit = dict(arcs=arcs, ann=ann, order=order, spans=spans, rframe=rframe, rann=rann, args=args, roots=roots,
               rank=rank)
    return prof, wit


def select(g, prof, wit, trial_a=0.00065):
    h, v = prof["h"], prof["v"]
    args, roots, order, ann, arcs = wit["args"], wit["roots"], wit["order"], wit["ann"], wit["arcs"]
    incoming = set(arcs.values())
    uses = [[] for _ in args]
    for x in order:
        if args[x]:
            for j, y in enumerate(args[x]):
                uses[y].append(2 * x + j)
    for j, r in enumerate(roots):
        uses[r["node"]].append(ROOTBIT | j)
    assign, sources, ops, coef, R = {}, {}, [], [], 0
    for x in order:
        if args[x]:
            aa, bb = args[x]
            dest, ctrl = assign[2 * x], assign[2 * x + 1]
            swapped = False
            if x in arcs:                      # the node's value is carried in the slot of one operand
                u = arcs[x]
                val = roots[u & (ROOTBIT - 1)]["node"] if u >> 31 else args[u // 2][u & 1]
                if val == aa:
                    dest, ctrl, swapped = ctrl, dest, True
                else:
                    assert val == bb
                assert u not in assign
                assign[u] = ctrl
            ops.append((dest, ctrl, x))
            sg = g["signs"][x - 1]
            coef.append((sg, 1) if swapped else (1, sg))
        else:
            dest = R
            R += 1
            sources[x] = dest
        free = [u for u in uses[x] if u not in incoming]
        assert free
        for j, u in enumerate(free):
            assert u not in assign
            if j == 0:
                assign[u] = dest
            else:                              # fan-out copy into a fresh role
                assign[u] = R
                ops.append((R, dest, x))
                coef.append((1, 1))
                R += 1
    assert R == prof["R"], (R, prof["R"])
    rootroles = [assign[ROOTBIT | j] for j in range(len(roots))]
    prev, pred, first = [-1] * R, [], [None] * R
    for i, (a, b, x) in enumerate(ops):
        pred.append((prev[a], prev[b]))
        prev[a] = prev[b] = i
        if first[a] is None:
            first[a] = x
        if first[b] is None:
            first[b] = x
    st = [prev[s] for r, s in zip(roots, rootroles) if r["kind"] == "center" and prev[s] >= 0]
    phase = set()
    while st:
        i = st.pop()
        if i not in phase:
            phase.add(i)
            st.extend(j for j in pred[i] if j >= 0)
    touched = set(sources.values())
    for i in phase:
        touched.update(ops[i][:2])
    co = [0] * R
    for r, s in zip(roots, rootroles):
        co[s] |= sum(1 << t for t in r["targets"])
    for a, b, x in reversed(ops):
        co[b] |= co[a]
    inputs, limit = g["inputs"], [None] * v
    for r in roots:
        if r["kind"] == "center":
            continue
        A = rref(inputs[t] for t in r["targets"])
        for t in r["targets"]:
            if limit[t] is None:
                limit[t] = A
    for t in range(v):
        if limit[t] is None:
            limit[t] = (inputs[t],)
    H = Counter(prof["H"])
    selected, gauges = [], Counter()
    cand = sorted((s for s in range(R) if s not in touched and first[s] is not None),
                  key=lambda s: (len(ann[first[s]]), bin(co[s]).count("1"), s))

    def excess(t):
        return t * math.expm1(trial_a * math.log(3 * h / t)) if t else 0.

    for s in cand:
        A, targets, bits = tuple(ann[first[s]]), [], co[s]
        while bits:
            low = bits & -bits
            t = low.bit_length() - 1
            bits ^= low
            targets.append(t)
            A = rref(A + tuple(limit[t]))
            if len(A) == h:
                break
        if len(A) == h:
            continue
        r, d = h - len(ann[first[s]]), h - len(A)
        delta = 3 * (excess(r - d) - excess(r)) + excess(3 * d)
        delta += 3 * math.fsum(excess(d) + excess(h - len(limit[t]) - d) - excess(h - len(limit[t])) for t in targets)
        if delta >= -1e-12:
            continue
        H[r] -= 1
        H[r - d] += 1
        gauges[d] += 1
        assert H[r] >= 0
        for t in targets:
            limit[t] = A
        selected.append(dict(role=s, annihilator=A, rank=d, targets=targets))
    word = dict(selected=selected, phase1=sorted(phase), sources=sources, rootroles=rootroles, ops=ops, coef=coef,
                R=R, gauges=gauges, H=H)
    return word
