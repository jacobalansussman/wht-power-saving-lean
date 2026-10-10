#!/usr/bin/env python3
"""comp.py -- compile a hierarchy of partial sums ("symbols") into an e8.Design.

Input (a plain object with these fields, e.g. gram.Grammar or a Hier below):
   content[u] = {T: coef}       what symbol u is (u < 120: the single source u)
   frame[u]                     span of the labels of its sources (echelon tuple)
   kids[u]    = None | [(child, coef), ...]     u = sum coef * child   (binary (u, v, r) of gram.py also accepted)
   need[S]    = {u: coef}       y_S must lose sum coef * u
Register allocation: every use of a symbol (by a parent rule at the parent's frame, or by a target at hyper(S)) needs a
register that holds the symbol and stands at that frame.  Uses at nested frames share one register (a chain).  A rule
is formed in place on the register of one child whose chain ends there; otherwise in a fresh helper.
Phase A: one carrier per chain of every source, loaded at the source's line; then the 8 totals at the full frame.
Phase B: the rules, by growing frame dimension; then the deliveries.  With ylow=True a target may take a symbol at the
frame where a register holding it already stands (one such chain of frames per target)."""
import sys
sys.dont_write_bytecode = True
from lib import *
from collections import defaultdict, Counter

def kids_of(h, u):
    k = h.kids[u]
    if k is None: return None
    if isinstance(k, tuple) and len(k) == 3 and not isinstance(k[0], tuple):
        return [(k[0], Fr(1)), (k[1], Fr(k[2]))]
    return [(c, Fr(x)) for c, x in k]

def compile_design(h, name="", ylow=False, stats=None, pairs=True, drop=frozenset()):
    n = len(h.content)
    HYP = [e8.hyper(S) for S in range(V)]
    # ---- which symbols are needed
    needed = set()
    todo = [u for S in range(V) for u in h.need[S]]
    while todo:
        u = todo.pop()
        if u in needed: continue
        needed.add(u)
        k = kids_of(h, u)
        if k: todo.extend(c for c, _ in k)
    order = sorted((u for u in needed if h.kids[u] is not None), key=lambda u: (len(h.frame[u]), u))
    # ---- uses: (dim, frame, kind, who)
    uses = defaultdict(list)
    for w in order:
        for c, _ in kids_of(h, w):
            uses[c].append((len(h.frame[w]), h.frame[w], 0, w))
    low = defaultdict(list)            # u -> [(S, coef)]: targets that take u at frame[u], where it is formed
    lowset = set()
    if ylow:
        for S in range(V):
            it = sorted(h.need[S], key=lambda u: (len(h.frame[u]), u))
            best = {}
            for j, u in enumerate(it):
                b = [u]
                for i in range(j):
                    v = it[i]
                    if (h.frame[v] == h.frame[u] or (len(h.frame[v]) < len(h.frame[u]) and e8.inside(h.frame[v], h.frame[u]))) \
                            and len(best[v]) + 1 > len(b):
                        b = best[v] + [u]
                best[u] = b
            ch = max(best.values(), key=len) if best else []
            if ylow != "all" and len(ch) < int(ylow):
                ch = []
            for u in ch:
                if len(h.frame[u]) < 8 and (S, u) not in drop:
                    low[u].append((S, h.need[S][u])); lowset.add((S, u))
    for S in range(V):
        for u in h.need[S]:
            if (S, u) not in lowset:
                uses[u].append((8, HYP[S], 1, S))
    # ---- chains (greedy: an event joins the chain whose last frame is the largest one inside its frame)
    chain_of = {}                      # (u, kind, who) -> chain index
    chains = {}                        # u -> list of chains, each a list of events
    for u in needed:
        ev = sorted(uses[u], key=lambda e: (e[0], e[2], e[3]))
        cl = []
        for e in ev:
            best = None
            for i, ch in enumerate(cl):
                last = ch[-1]
                if last[2] == 1 and last[1] != e[1]:
                    continue                           # after a target read at hyper(S) nothing else
                if last[1] == e[1] or (last[0] < e[0] and e8.inside(last[1], e[1])):
                    if best is None or cl[best][-1][0] < last[0]:
                        best = i
            if best is None:
                cl.append([e]); best = len(cl) - 1
            else:
                cl[best].append(e)
            chain_of[(u, e[2], e[3])] = best
        chains[u] = cl
    d = e8.Design(name=name or "comp.py")
    reg = {}                           # (u, chain) -> (register, scale): register holds scale * content[u]
    # ---- phase A: carriers
    for T in range(V):
        if T in needed:
            if not chains[T]: chains[T].append([])
            qs = d.helpers(len(chains[T]))
            d.out("A", e8.line(T), d.X(T), [(q, 1) for q in qs])
            for i, q in enumerate(qs): reg[(T, i)] = (q, Fr(1))
    cen = e8.totals8()
    tot = d.helpers(len(cen))
    for k, (r, g) in enumerate(cen):
        d.inn("A", e8.FULL, tot[k], [(d.X(T), x) for T, x in sorted(g.items()) if x])
        d.retain(tot[k])
    d.scatter([[(k, r[S]) for k, (r, g) in enumerate(cen) if r.get(S, 0)] for S in range(V)])
    # ---- phase B: rules
    inplace = 0
    twin = getattr(h, "twin", {}) if pairs else {}
    done_pair = set(); npair = 0
    for T in range(V):
        if T in low:
            q, sc = reg[(T, 0)]
            d.out("B", e8.line(T), q, [(d.Y(S), -coef / sc) for S, coef in low[T]])
    for w in order:
        F = h.frame[w]
        ks = kids_of(h, w)
        if w in done_pair:
            continue
        w2 = twin.get(w)
        if w2 is not None and w2 in needed and len(ks) == 2:
            # butterfly: u + r v and u + r2 v, both formed in place on the registers of u and v
            (u, xu), (v, xv) = ks
            (u2, yu), (v2, yv) = kids_of(h, w2)
            if (u2, v2) == (u, v) and xu == 1 and yu == 1 and xv != yv:
                iu, iv = chain_of[(u, 0, w)], chain_of[(v, 0, w)]
                lu = [e[3] for e in chains[u][iu] if e[2] == 0 and e[1] == F]; lv = [e[3] for e in chains[v][iv] if e[2] == 0 and e[1] == F]
                tail_u = chains[u][iu][-len(lu):] if lu else []
                tail_v = chains[v][iv][-len(lv):] if lv else []
                okp = (chain_of[(u, 0, w2)] == iu and chain_of[(v, 0, w2)] == iv and sorted(lu) == sorted([w, w2]) and sorted(lv) == sorted([w, w2])
                       and len(tail_u) == 2 and len(tail_v) == 2 and all(e[2] == 0 for e in tail_u + tail_v)
                       and chains[u][iu][-1][1] == F and chains[v][iv][-1][1] == F)
                if okp:
                    qu, su = reg[(u, iu)]; qv, sv = reg[(v, iv)]
                    # THE IN-PLACE PAIR (a + b and a - b formed on the two arrays that held a and b) is an outside device:
                    # icekylinx, hub #184 (CrocSwap/integer-mult-bounds), carried in ikeboy's #191 and #193.
                    # the first version of this rule was abs(su) == 1 and abs(sv) == 1.  One of the two helpers may hold halves:
                    # then both sums come out in halves (never quarters).  Two helpers in halves stay excluded.
                    if (abs(su), abs(sv)) in ((1, 1), (Fr(1, 2), 1), (1, Fr(1, 2))) and yv == -xv:
                        reg.pop((u, iu)); reg.pop((v, iv))
                        if abs(sv) == 1:
                            d.inn("B", F, qu, [(qv, xv * su / sv)])                 # qu = su (u + xv v) = su w
                            cc = sv / (su * (yv - xv))
                            d.inn("B", F, qv, [(qu, cc)])                          # qv = sv/(yv - xv) (u + yv v)
                            outs = ((w, qu, su), (w2, qv, sv / (yv - xv)))
                        else:
                            c1 = sv / (xv * su)
                            d.inn("B", F, qv, [(qu, c1)])                          # qv = (sv/xv) (u + xv v) = (sv/xv) w
                            c2 = su * yv / (2 * sv)
                            d.inn("B", F, qu, [(qv, c2)])                          # qu = (su/2) (u + yv v) = (su/2) w2
                            outs = ((w, qv, sv / xv), (w2, qu, su / 2))
                        for ww, q, sc in outs:
                            nch = len(chains[ww])
                            reg[(ww, 0)] = (q, sc)
                            if nch > 1:
                                qs = d.helpers(nch - 1)
                                d.out("B", F, q, [(qq, 1 / sc) for qq in qs])
                                for i, qq in enumerate(qs): reg[(ww, i + 1)] = (qq, Fr(1))
                            if ww in low:
                                d.out("B", F, q, [(d.Y(S), -coef / sc) for S, coef in low[ww]])
                        done_pair.add(w2); npair += 1; inplace += 2
                        continue
        ch = []
        for c, x in ks:
            i = chain_of[(c, 0, w)]
            last = chains[c][i][-1]
            ch.append((c, x, i, last[2] == 0 and last[3] == w))
        host = next((t for t in ch if t[3]), None)
        if host is not None:
            c0, x0, i0, _ = host
            q0, s0 = reg.pop((c0, i0))
            # register holds s0 * c0 ; want lambda * w with lambda = s0 / x0 :  add (s0/x0) * x * child
            lam = s0 / x0
            srcs = []
            for c, x, i, _ in ch:
                if (c, i) == (c0, i0): continue
                q, s = reg[(c, i)]
                srcs.append((q, lam * x / s))
            d.inn("B", F, q0, srcs)
            wreg, wsc = q0, lam
            inplace += 1
        else:
            wreg = d.helper(); wsc = Fr(1)
            d.inn("B", F, wreg, [(reg[(c, i)][0], x / reg[(c, i)][1]) for c, x, i, _ in ch])
        nch = len(chains[w])
        reg[(w, 0)] = (wreg, wsc)
        if nch > 1:
            qs = d.helpers(nch - 1)
            d.out("B", F, wreg, [(q, 1 / wsc) for q in qs])
            for i, q in enumerate(qs): reg[(w, i + 1)] = (q, Fr(1))
        if w in low:
            d.out("B", F, wreg, [(d.Y(S), -coef / wsc) for S, coef in low[w]])
    # ---- deliveries
    for S in range(V):
        srcs = []
        for u, coef in sorted(h.need[S].items()):
            if (S, u) in lowset: continue
            q, s = reg[(u, chain_of[(u, 1, S)])]
            srcs.append((q, -coef / s))
        if srcs: d.inn("B", HYP[S], d.Y(S), srcs)
    if stats is not None:
        stats.update(R=d.R, rules=len(order), inplace=inplace, carriers=sum(len(chains[T]) for T in range(V) if T in needed),
                     chains=sum(len(c) for c in chains.values()), lowreads=len(lowset), butterflies=npair, lowset=sorted(lowset))
    return d
