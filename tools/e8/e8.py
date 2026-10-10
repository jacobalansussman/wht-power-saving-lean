#!/usr/bin/env python3
"""e8.py -- shared module of the E8 search programs.  Python 3.9 or later, standard library only.

What is here
  (a) the E8 family: NAMES, PORTS, VEC, knows(S, T), BMAT (fitting matrix), totals8(), INCIDENCES
      and subspace arithmetic over F_2 (echelon, perp, join, meet, inside, line, hyper, FULL, ZERO);
  (b) Design + build_cert(design) + write_cert(design, path): a design is exactly what the certificate format
      leaves free (number of helpers, the two gate lists with their frames and coefficients, the retained totals,
      the scatter).  Everything else in the file is forced or derived and is filled in here.
      Convenience builders: partial_sum, total, carrier, carriers, collector, deliver, naive_design;
  (c) quick(cert)      IN-PROCESS: the reference checker gx.check1 and the mirror gxcore.mirror of ../gx, called as
                       functions, plus the predicted figure.  This is what the search loops use.  It is not the
                       verdict: the verdict is  python3 -B tools/gx/refcheck.py <file>  and  tools/gx/gxdry.py <file>
                       (run.py calls both on the final file);
  (e) predict(R, block_profile, cst): the whole-block figure of a hypothetical circuit, by the code of
      tools/gx/gxunit.py (unit, rates) and tools/gen/foldlib.py (best_units, plan_block, plan_rank).
"""
import sys
sys.dont_write_bytecode = True
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import os
import re
import time
from collections import Counter
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))            # tools/e8
TOOLS = os.path.dirname(HERE)                                # tools
GX, GEN = TOOLS + "/gx", TOOLS + "/gen"                      # the published checkers and generators, used unchanged

# =====================================================================================================================
# (a) the E8 family
# =====================================================================================================================
H = 9                      # label width: labels are bit vectors of length 9 (integers 1 .. 511; bit i = point i)
NPOINTS = 9
FULLMASK = (1 << H) - 1


def par(x):
    """parity of the number of one-bits"""
    return bin(x).count("1") & 1


def family(n=NPOINTS):
    """names, labels (ints), lattice vectors (a_1..a_n ; t).
    Order: first the C(n,3) triples in the order of itertools.combinations(range(n), 3) (lexicographic),
    then the C(n,2) pairs in the order of itertools.combinations(range(n), 2)."""
    names, lab, vec = [], [], []
    full = (1 << n) - 1
    for T in itertools.combinations(range(n), 3):
        names.append(("T",) + T)
        lab.append(sum(1 << i for i in T))
        vec.append([1 if i in T else 0 for i in range(n)] + [1])
    for i, j in itertools.combinations(range(n), 2):
        names.append(("P", i, j))
        lab.append(full ^ (1 << i) ^ (1 << j))
        a = [0] * (n + 1)
        a[i] = 1
        a[j] = -1
        vec.append(a)
    return names, lab, vec


NAMES, PORTS, VEC = family(NPOINTS)
V = len(PORTS)             # 120 pairs: indices 0..83 = triples (weight-3 labels), 84..119 = pairs (weight-7 labels)
NTRIPLES = 84
PORTSET = frozenset(PORTS)
INDEX = {u: i for i, u in enumerate(PORTS)}


def knows(S, T):
    """S knows T: S != T and the labels have even overlap (orthogonal over F_2).  Symmetric.
    Then line(T) lies inside hyper(S), so the content of x_T can reach y_S by climbing only."""
    return S != T and not par(PORTS[S] & PORTS[T])


def gram():
    """table of products <w, w'> = a.a' - t t' of the lattice vectors"""
    n = NPOINTS
    return [[sum(a[k] * b[k] for k in range(n)) - a[n] * b[n] for b in VEC] for a in VEC]


GRAM = gram()
BMAT = [[Fr(g, 2) for g in row] for row in GRAM]     # fitting matrix: B[S][S] = 1, B[S][T] = +-1/2 iff S knows T, else 0
KNOWN = [[T for T in range(V) if knows(S, T)] for S in range(V)]         # 56 per target
INCIDENCES = [(S, T) for S in range(V) for T in KNOWN[S] if BMAT[S][T]]   # 6,720 ordered pairs (target, source)


def totals8():
    """the 8 known totals of E8, as [(r_k, g_k)], k = 0..7:
       g_k[T] = a_k(T) - a_8(T)          what total k holds: total_k = sum_T g_k[T] x_T
       r_k[S] = (a_k(S) - t(S)/3) / 2    the scatter coefficient: y_S += r_k[S] * total_k
    and sum_k r_k[S] g_k[T] = BMAT[S][T].  With all 8 totals standing at the full frame, cst = 8 * 9 = 72."""
    n = NPOINTS
    cen = []
    for i in range(n - 1):
        g = {T: Fr(VEC[T][i] - VEC[T][n - 1]) for T in range(V) if VEC[T][i] - VEC[T][n - 1]}
        r = {S: (Fr(VEC[S][i]) - Fr(VEC[S][n], 3)) / 2 for S in range(V)}
        cen.append(({S: x for S, x in r.items() if x}, g))
    return cen


CST8 = 8 * H               # 72


# ---- subspaces of F_2^h.  A frame is a tuple of ints: reduced echelon basis, leading bits descending ----------------
def echelon(vs):
    """reduced echelon basis of the span.  This is the normal form the checkers require
    (gx.py:30-34: no zero vector, strictly descending leading bits, each leading bit in one basis vector only)."""
    basis = []
    for u in vs:
        for b in basis:
            if u >> (b.bit_length() - 1) & 1:
                u ^= b
        if u:
            lb = u.bit_length() - 1
            basis = [b ^ u if b >> lb & 1 else b for b in basis]
            basis.append(u)
            basis.sort(reverse=True)
    return tuple(basis)


def reduce(u, basis):
    for b in basis:
        if u >> (b.bit_length() - 1) & 1:
            u ^= b
    return u


def inside(a, b):
    """span(a) is contained in span(b); b must be in echelon form"""
    return all(reduce(u, b) == 0 for u in a)


def perp(basis, h=H):
    """{x in F_2^h : x.b = 0 for all b in basis}, in echelon form"""
    Bs = echelon(basis)
    piv = [b.bit_length() - 1 for b in Bs]
    out = []
    for k in range(h):
        if k in piv:
            continue
        x = 1 << k
        for b, p in zip(Bs, piv):
            if b >> k & 1:
                x |= 1 << p
        out.append(x)
    return echelon(out)


def join(*frames):
    """smallest subspace containing all the given ones"""
    return echelon([u for f in frames for u in f])


def meet(f1, f2, h=H):
    """intersection"""
    return perp(tuple(perp(f1, h)) + tuple(perp(f2, h)), h)


ZERO = ()
FULL = tuple(1 << i for i in range(H - 1, -1, -1))


def line(T):
    """the frame where x_T starts: the line spanned by its label"""
    return (PORTS[T],)


def hyper(S):
    """the frame where y_S must end: all vectors with even overlap with the label of S (dimension h - 1)"""
    return perp((PORTS[S],), H)


# =====================================================================================================================
# (b) designs and certificates
# =====================================================================================================================
class DesignError(Exception):
    pass


class NoSlot(RuntimeError):
    """not used by the published programs"""


def _co(c):
    """coefficient -> [numerator, denominator], denominator > 0, reduced"""
    if isinstance(c, Fr):
        f = c
    elif isinstance(c, int):
        f = Fr(c)
    elif isinstance(c, (tuple, list)) and len(c) == 2:
        f = Fr(int(c[0]), int(c[1]))
    else:
        raise DesignError("coefficient %r: give an int, a Fraction or a pair (num, den)" % (c,))
    return [f.numerator, f.denominator]


def _items(x):
    """{register: coef}, or a list of (register, coef), or a list of [register, num, den] as in the certificate"""
    if isinstance(x, dict):
        return list(x.items())
    return [(e[0], (e[1], e[2])) if len(e) == 3 else (e[0], e[1]) for e in x]


class Design:
    """Everything the gcert/1 format leaves free, and nothing else.

    registers   X(T) = T, Y(S) = v + S, helper q = 2 v + q  (the numbering of the certificate; gx.py:148).
    frames      any iterable of ints spanning the subspace (it is brought to echelon form here), or FULL / ZERO.
    gates       ("out", frame, src, [(tgt, coef)...], extra)   every tgt += coef * src (src read once, at the start)
                ("in",  frame, tgt, [(src, coef)...])          tgt += sum coef * src, one after the other
                Every register a gate names (extras too) first climbs to the gate's frame.  A gate with no adds
                ("out" with an empty target list) is a pure STOP: see stop().
    phase       "A" = before the scatter, "B" = after it.  The order inside each list is the order of execution.
    ret         the helpers whose contents, at the moment of the scatter, are total 0, 1, 2, ...  The frame of
                total k is wherever that helper stands then; its scratch copy costs dim(frame) moves (this is cst).
    scat        table: one row per target S, [(k, coef)...]:  y_S += coef * total_k  (all y at frame 0);
                or the star rule (inside, outside), which needs exactly h totals (gx.py:116-117, 171).
    Forced, not free: x_T starts on line(T) and ends at FULL; y_S starts at ZERO and ends at hyper(S); every helper
    starts at ZERO and ends at FULL (gx.py:142-147); blocks, N and cst are recomputed by the checkers.
    """

    def __init__(self, ports=None, h=H, name=""):
        self.h = h
        self.ports = list(PORTS if ports is None else ports)
        self.v = len(self.ports)
        self.R = 0
        self.A, self.B = [], []
        self.ret = []
        self.scat = None
        self.name = name
        self._fc = {}

    # registers
    def X(self, T):
        return T

    def Y(self, S):
        return self.v + S

    def helper(self):
        """a new helper register (starts at frame 0 with no x-content)"""
        self.R += 1
        return 2 * self.v + self.R - 1

    def helpers(self, n):
        return [self.helper() for _ in range(n)]

    # frames
    def frame(self, f):
        key = tuple(f)
        e = self._fc.get(key)
        if e is None:
            e = self._fc[key] = echelon(key)
        return e

    def _ph(self, phase):
        if phase not in ("A", "B"):
            raise DesignError("phase must be 'A' or 'B'")
        return self.A if phase == "A" else self.B

    # gates
    def out(self, phase, frame, src, targets, extra=()):
        """fan-out: every target += coef * src.  targets: dict {register: coef} or list of (register, coef)."""
        g = ["out", self.frame(frame), int(src), [[int(t)] + _co(c) for t, c in _items(targets)], [int(r) for r in extra]]
        self._ph(phase).append(g)
        return g

    def inn(self, phase, frame, tgt, sources):
        """fan-in: tgt += sum coef * src.  sources: dict {register: coef} or list of (register, coef)."""
        g = ["in", self.frame(frame), int(tgt), [[int(s)] + _co(c) for s, c in _items(sources)]]
        self._ph(phase).append(g)
        return g

    def stop(self, phase, frame, regs):
        """a gate without adds: the registers climb to `frame` and stop there (one block ends for each)"""
        regs = [int(r) for r in regs]
        g = ["out", self.frame(frame), regs[0], [], regs[1:]]
        self._ph(phase).append(g)
        return g

    # scatter
    def retain(self, slot):
        """the helper `slot` holds total k at the scatter; returns k"""
        self.ret.append(int(slot))
        return len(self.ret) - 1

    def scatter(self, table):
        """table[S] = {k: coef} or [(k, coef)...]"""
        self.scat = ("table", [[[int(k)] + _co(c) for k, c in _items(row)] for row in table])

    def scatter_star(self, inside, outside):
        self.scat = ("star", _co(inside), _co(outside))


# ---- convenience builders (thin: each adds one or two gates) -----------------------------------------------------
def partial_sum(d, coeffs, frame=FULL, phase="A"):
    """a new helper that holds sum_T coeffs[T] * x_T, formed by ONE fan-in gate at `frame`.
    `frame` must contain line(T) for every T named; every x_T named climbs to `frame` (and can never come back)."""
    q = d.helper()
    d.inn(phase, frame, q, [(d.X(T), c) for T, c in _items(coeffs)])
    return q


def total(d, coeffs, frame=FULL):
    """partial_sum in phase A, retained as the next total.  Returns (helper, k)."""
    q = partial_sum(d, coeffs, frame, "A")
    return q, d.retain(q)


def carriers(d, T, n, frame=None, phase="A"):
    """n new single-source carriers of x_T, loaded by ONE fan-out gate at `frame` (default: line(T), where x_T
    starts, so that x_T does not move).  Returns the list of helpers; each holds 1 * x_T."""
    qs = d.helpers(n)
    d.out(phase, line(T) if frame is None else frame, d.X(T), [(q, 1) for q in qs])
    return qs


def carrier(d, T, frame=None, phase="A"):
    return carriers(d, T, 1, frame, phase)[0]


def collector(d, sources, frame, phase="B"):
    """a new helper that receives sum coef * register over `sources` ({register: coef}; x registers or helpers)
    by ONE fan-in gate at `frame`.  All the sources climb to `frame`."""
    q = d.helper()
    d.inn(phase, frame, q, sources)
    return q


def deliver(d, S, sources, frame=None):
    """phase B: y_S += sum coef * register over `sources`, at `frame` (default hyper(S), the last frame of y_S).
    `frame` must lie inside hyper(S); x registers and helpers may be read, y registers only under rule E4."""
    return d.inn("B", hyper(S) if frame is None else frame, d.Y(S), sources)


def naive_design():
    """the WASTEFUL scheme, written through this module: one carrier per incidence.
    Returns (design, slot, tot): slot[S, T] = the carrier of source T for target S; tot[k] = helper of total k."""
    d = Design(name="e8.naive_design: the wasteful scheme (one carrier per incidence)")
    v = d.v
    cen = totals8()
    slot = {}
    for S in range(v):
        for T in KNOWN[S]:
            if BMAT[S][T]:
                slot[S, T] = d.helper()
    tot = d.helpers(len(cen))
    for T in range(v):
        tg = [(slot[S, T], 1) for S in range(v) if (S, T) in slot]
        if tg:
            d.out("A", line(T), d.X(T), tg)
    for k, (r, g) in enumerate(cen):
        d.inn("A", FULL, tot[k], [(d.X(T), x) for T, x in sorted(g.items()) if x])
        d.retain(tot[k])
    d.scatter([[(k, r[S]) for k, (r, g) in enumerate(cen) if r.get(S, 0)] for S in range(v)])
    for S in range(v):
        deliver(d, S, [(slot[S, T], -BMAT[S][T]) for T in KNOWN[S] if (S, T) in slot])
    return d, slot, tot


# ---- design -> certificate -----------------------------------------------------------------------------------------
def gate_regs(g):
    """registers named by a certificate gate, in the order they move (Raw.lean Gate.regs; gx.py:120-126)"""
    if g[0] == "out":
        return [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    return [g[2]] + [s for s, _, _ in g[3]]


def _as_design(design):
    if isinstance(design, Design):
        return design
    d = Design(ports=design.get("ports"), h=design.get("h", H), name=design.get("name", ""))
    d.R = int(design["R"])
    for ph in "AB":
        for g in design[ph]:
            if g[0] == "out":
                d.out(ph, g[1], g[2], g[3], g[4] if len(g) > 4 else ())
            elif g[0] == "in":
                d.inn(ph, g[1], g[2], g[3])
            else:
                raise DesignError("unknown gate kind %r" % (g[0],))
    for s in design["ret"]:
        d.retain(s)
    sc = design["scat"]
    if isinstance(sc, dict) and "inside" in sc:
        d.scatter_star(sc["inside"], sc["outside"])
    else:
        d.scatter(sc["table"] if isinstance(sc, dict) else sc)
    return d


def build_cert(design, strict=True):
    """Design (or a dict with keys h, ports, R, A, B, ret, scat in the conventions of Design) -> gcert/1 dict.
    Fills in the forced and derived fields: frame table, start, final, ret frames, blocks, N, cst, ext = [].
    strict=True raises DesignError with a readable message when a register would have to move to a frame that does
    not strictly contain its present one (the checkers' rule E1); strict=False builds the file anyway."""
    d = _as_design(design)
    h, ports, v, R = d.h, d.ports, d.v, d.R
    if d.scat is None:
        raise DesignError("no scatter given")
    full = tuple(1 << i for i in range(h - 1, -1, -1))
    frames, index = [(), full], {(): 0, full: 1}

    def fid(b):
        i = index.get(b)
        if i is None:
            i = index[b] = len(frames)
            frames.append(b)
        return i
    lines = [fid((u,)) for u in ports]
    hyp = [fid(perp((u,), h)) for u in ports]
    conv = lambda g: ["out", fid(g[1]), g[2], g[3], g[4]] if g[0] == "out" else ["in", fid(g[1]), g[2], g[3]]
    A, Bg = [conv(g) for g in d.A], [conv(g) for g in d.B]
    nreg = 2 * v + R
    start = lines + [0] * v + [0] * R
    final = [1] * v + hyp + [1] * R
    dim = [len(f) for f in frames]
    cur = list(start)
    Hh = {"x": {}, "y": {}, "s": {}, "c": {}}
    cls = lambda r: "x" if r < v else "y" if r < 2 * v else "s"
    ok = {}

    def climb(r, f, where):
        o = cur[r]
        if o != f:
            if strict:
                if not 0 <= r < nreg:
                    raise DesignError("%s: register %d does not exist (2v + R = %d)" % (where, r, nreg))
                if (o, f) not in ok:
                    ok[o, f] = dim[o] < dim[f] and inside(frames[o], frames[f])
                if not ok[o, f]:
                    raise DesignError("%s: register %d (%s) stands at frame %s (dim %d) and cannot climb to %s (dim %d)"
                                      % (where, r, cls(r), list(frames[o]), dim[o], list(frames[f]), dim[f]))
            k = cls(r) if 0 <= r < nreg else "s"
            Hh[k][dim[f] - dim[o]] = Hh[k].get(dim[f] - dim[o], 0) + 1
            cur[r] = f
    for i, g in enumerate(A):
        for r in gate_regs(g):
            climb(r, g[1], "gate A%d" % i)
    ret = [[k, s, cur[s]] for k, s in enumerate(d.ret)]
    for k, s, f in ret:
        Hh["c"][dim[f]] = Hh["c"].get(dim[f], 0) + 1
    for i, g in enumerate(Bg):
        for r in gate_regs(g):
            climb(r, g[1], "gate B%d" % i)
    for r in range(nreg):
        climb(r, final[r], "final climb")
    blocks = {k: {str(r): n for r, n in sorted(w.items())} for k, w in Hh.items()}
    N = sum(r * n for w in Hh.values() for r, n in w.items())
    cst = sum(dim[f] for _, _, f in ret)
    if d.scat[0] == "table":
        scat = dict(table=d.scat[1])
    else:
        scat = dict(inside=d.scat[1], outside=d.scat[2])
    return dict(format="gcert/1", derived_from=d.name or "e8.py", p=0, h=h, v=v, R=R, cst=cst, N=N, ports=list(ports),
                frames=[list(f) for f in frames], start=start, final=final, ext=[], ret=ret, scat=scat, A=A, B=Bg,
                blocks=blocks)


def dump(c, path):
    """write a certificate dict; .gz files are written with a zero time stamp, so equal contents give equal sha256"""
    data = json.dumps(c, separators=(",", ":")).encode()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if path.endswith(".gz"):
        with open(path, "wb") as fh:
            with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0) as gz:
                gz.write(data)
    else:
        with open(path, "wb") as fh:
            fh.write(data)


def load(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        return json.load(fh)


def write_cert(design, path, strict=True):
    """design -> gcert/1 file at `path` (.json or .json.gz).  Returns the certificate dict."""
    c = build_cert(design, strict)
    dump(c, path)
    return c


def sha256(path):
    hsh = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            hsh.update(blk)
    return hsh.hexdigest()


# ---- reading a certificate back ------------------------------------------------------------------------------------
def paths(c):
    """for every register: (class, [frame ids it stops at, start included], [block ranks]) by plain replay"""
    v, R = c["v"], c["R"]
    dim = [len(f) for f in c["frames"]]
    cur = list(c["start"])
    stops = [[f] for f in cur]
    for g in c["A"] + c["B"]:
        for r in gate_regs(g):
            if cur[r] != g[1]:
                cur[r] = g[1]
                stops[r].append(g[1])
    for r, f in enumerate(c["final"]):
        if cur[r] != f:
            stops[r].append(f)
    cls = lambda r: "x" if r < v else "y" if r < 2 * v else "s"
    return [(cls(r), st, [dim[b] - dim[a] for a, b in zip(st, st[1:])]) for r, st in enumerate(stops)]


def profile(c):
    """{class: Counter of block patterns}: how the registers of each class cut their climb into blocks"""
    out = {"x": Counter(), "y": Counter(), "s": Counter()}
    for k, _, ranks in paths(c):
        out[k][tuple(ranks)] += 1
    return out


# =====================================================================================================================
# the published tools (../gx, ../gen)
# =====================================================================================================================
_MODS = {}


def _tool(name):
    """import a module of the published tools (gx, gxcore from tools/gx; foldlib from tools/gen) without writing anything"""
    if name not in _MODS:
        if name == "foldlib":
            sys.path.insert(0, GEN)
            try:
                import foldlib
            finally:
                sys.path.remove(GEN)
            _MODS[name] = foldlib
        else:
            spec = importlib.util.spec_from_file_location("e8tools_" + name, "%s/%s.py" % (GX, name))
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            _MODS[name] = m
    return _MODS[name]


# =====================================================================================================================
# (e) the figure
# =====================================================================================================================
def unit_blocks(inv, h, v, R, cst):
    """block list of one five-stage unit, copied from tools/gx/gxunit.py:22-36:
    5 x (all blocks of one invocation: x, y, helpers, scratch copies) + per pair 2 blocks each of rank
    h-1, 2h-2, 2h+2, 4 (tools/gen/foldlib.py:24, BANK['B2']).  m = 5 h, W = 4 v + R, D = W m - moves."""
    fl = _tool("foldlib")
    Hh = {}

    def add(r, n):
        if n:
            Hh[r] = Hh.get(r, 0) + n
    for r, n in inv.items():
        add(r, 5 * n)
    for cv, rk in fl.BANK['B2']:
        add(rk(h), cv * v)
    m, W = 5 * h, 4 * v + R
    D = W * m - sum(r * n for r, n in Hh.items())
    if D != 4 * v - 5 * cst:
        raise DesignError("block list does not add up: D = %d but 4v - 5 cst = %d (helper moves must total R h, "
                          "x and y moves v (h-1) each, copies cst)" % (D, 4 * v - 5 * cst))
    return dict(kind='B2', h=h, m=m, W=W, D=D, v=v, R=R, H=Hh, inv=dict(inv), cst=cst)


def rates(u, B=10, s=40):
    """copied from tools/gx/gxunit.py:39-47: the largest A (numerator over 10^B) at which the whole-block fact holds
    by 80-digit evaluation (foldlib.best_units), lowered until a rational certificate exists (foldlib.plan_block);
    the same for the per-rank fact."""
    fl = _tool("foldlib")
    A0, Af, Ar0, Arf = fl.best_units(u, B, s)
    A, Ar = A0, Ar0
    while fl.plan_block(u, A, B, s) is None:
        A -= 1
    while fl.plan_rank(u, Ar, B, s) is None:
        Ar -= 1
    return dict(A=A, A_eval=A0, A_fails=Af, A_rank=Ar, A_rank_eval=Ar0, A_rank_fails=Arf, B=B, s=s)


def _helper_hist(R, prof, h):
    """helper blocks of ONE invocation as {rank: count}"""
    if isinstance(prof, str):
        key = prof.replace(" ", "")
        if key in ("three", "1,h-2,1"):
            prof = (1, h - 2, 1)
        elif key == "unit":
            prof = (1,) * h
        elif key == "one":
            prof = (h,)
        elif key == "wasteful":
            if R < 8:
                raise DesignError("profile 'wasteful' needs R >= 8")
            prof = {(1, h - 2, 1): R - 8, (h,): 8}
        else:
            raise DesignError("unknown profile %r" % prof)
    if isinstance(prof, (tuple, list)):
        prof = {tuple(prof): R}
    hist = {}
    if prof and all(isinstance(k, tuple) for k in prof):          # {pattern: number of helpers}
        if sum(prof.values()) != R:
            raise DesignError("patterns cover %d helpers, R = %d" % (sum(prof.values()), R))
        for pat, n in prof.items():
            if sum(pat) != h or min(pat) < 1:
                raise DesignError("pattern %r: the block ranks of a helper must be positive and sum to h = %d" % (pat, h))
            for r in pat:
                hist[r] = hist.get(r, 0) + n
    else:                                                         # {rank: count}
        hist = {int(r): int(n) for r, n in prof.items() if n}
        if sum(r * n for r, n in hist.items()) != R * h:
            raise DesignError("helper blocks hold %d moves, must be R h = %d" % (sum(r * n for r, n in hist.items()), R * h))
    return hist


def predict(R, block_profile, cst=CST8, v=V, h=H, x_blocks=None, y_blocks=None, copy_ranks=None, B=10, s=40):
    """the figure gxdry.py would print for a circuit with R helpers.

    block_profile  how the helpers cut their h moves into blocks in ONE invocation:
                   'three' (= '1,h-2,1'), 'unit' (h blocks of rank 1), 'one' (one block of rank h),
                   'wasteful' (8 helpers of one block h, the others 1, h-2, 1);
                   a pattern such as (1, 7, 1) for every helper;  {pattern: number of helpers};
                   or {rank: count}, the histogram of all helper blocks (class "s" of the certificate).
    cst            copy moves per invocation;  copy_ranks = the frame dimensions of the totals (default: cst / h
                   totals at the full frame).
    x_blocks, y_blocks   {rank: count} of the x and y registers (default: one block of rank h - 1 each).
    Returns whole_block (the score), per_rank, m, W, D, the unit block list, and two float diagnostics:
    root (1 - z at M(z) = W, times 10^B) and first_order (D / sum of r ln(m/r) over all blocks, times 10^B)."""
    if copy_ranks is None:
        if cst % h:
            raise DesignError("cst = %d is not a multiple of h: give copy_ranks" % cst)
        copy_ranks = [h] * (cst // h)
    if sum(copy_ranks) != cst:
        raise DesignError("copy_ranks sum to %d, cst = %d" % (sum(copy_ranks), cst))
    inv = {}

    def add(hist):
        for r, n in hist.items():
            if n:
                inv[r] = inv.get(r, 0) + n
    add(x_blocks if x_blocks is not None else {h - 1: v})
    add(y_blocks if y_blocks is not None else {h - 1: v})
    add(_helper_hist(R, block_profile, h))
    add(Counter(copy_ranks))
    return figure(inv, h, v, R, cst, B, s)


def figure(inv, h, v, R, cst, B=10, s=40):
    """the figures of the five-stage unit from the block histogram of one invocation ({rank: count}, all classes)"""
    u = unit_blocks(inv, h, v, R, cst)
    r = rates(u, B, s)
    m, W, D = u["m"], u["W"], u["D"]
    C = sum(n * rk * math.log(m / rk) for rk, n in u["H"].items())
    f = lambda z: sum(n * (rk / m) ** z for rk, n in u["H"].items()) - W
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return dict(whole_block=r["A"], whole_block_eval=r["A_eval"], whole_block_fails=r["A_fails"], per_rank=r["A_rank"],
                m=m, W=W, D=D, R=R, cst=cst, unit=dict(sorted(u["H"].items())), inv=dict(sorted(inv.items())),
                root=(1 - (lo + hi) / 2) * 10 ** B if D > 0 else 0.0, first_order=D / C * 10 ** B,
                moves_weight=C, B=B, s=s)


def inv_hist(c):
    """block histogram of one invocation from the `blocks` field of a certificate"""
    inv = {}
    for k in "xysc":
        for r, n in c["blocks"][k].items():
            inv[int(r)] = inv.get(int(r), 0) + n
    return inv


# =====================================================================================================================
# (c) checking
# =====================================================================================================================
def quick(cert, rate=True):
    """IN-PROCESS check of a certificate (dict or path) by the two published checkers, called as functions:
    gx.check1 (what refcheck.py runs) and gxcore.mirror (what gxdry.py runs).  Not the official verdict: run
    refcheck.py and gxdry.py on the file for that.  Returns accepted_ref, accepted_mirror, the two messages, h v R N cst, and (if both accept and rate)
    the predicted whole_block / per_rank / m / W / D."""
    c = load(cert) if isinstance(cert, str) else json.loads(json.dumps(cert))
    gx, gc = _tool("gx"), _tool("gxcore")
    out = dict(accepted_ref=False, accepted_mirror=False, ref_msg="", mirror_msg="", h=c.get("h"), v=c.get("v"),
               R=c.get("R"), N=c.get("N"), cst=c.get("cst"), whole_block=None, per_rank=None)
    st = {}
    try:
        gx.check1(c, scalar=True, stats=st)
        out["accepted_ref"] = True
        out["ref_stats"] = st
    except AssertionError as e:
        out["ref_msg"] = str(e)
    except Exception as e:                      # a malformed file can crash the reference checker: that is a refusal
        out["ref_msg"] = "crash: %r" % (e,)
    try:
        D = gc.normal(c)
        gu_adds = sum(len(g[3]) for g in D['A'] + D['B']) + (D['v'] * D['h'] if D['scat'][0] == 'star'
                                                              else sum(map(len, D['scat'][1])))
        M = gc.mirror(D, gmax=12000, batches=1 if gu_adds <= 30000 else 2)      # gxunit.py:8-19 (plan)
        out["accepted_mirror"] = True
        out["hist"] = M["hist"]
        out["H"] = M["H"]
        out["adds"] = gu_adds
        out["most_regs_one_gate"] = M["lab"]["most"]
    except gc.Reject as e:
        out["mirror_msg"] = str(e)
    except Exception as e:
        out["mirror_msg"] = "crash: %r" % (e,)
    if rate and out["accepted_ref"] and out["accepted_mirror"] and not c.get("ext"):
        f = figure(out["hist"], c["h"], c["v"], c["R"], c["cst"])
        out.update(whole_block=f["whole_block"], per_rank=f["per_rank"], m=f["m"], W=f["W"], D=f["D"], unit=f["unit"])
    return out
