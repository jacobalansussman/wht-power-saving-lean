# Recombiner check (wht44/ideas, 16 Oct 2026): which ENGINE-GAME positions at block 4 pass the f = 9 hook's
# termwise filter (rho_9(g) E inside N, E = N cap K^-1 N)?
#
# Why it matters.  px-t514 Lemma A (engine game, G finite): a unit exactly tight on a certificate Y at one f is tight on
# Y (x) pi for every irreducible pi; in the regular representation this is Theorem T (X7): EVERY TERM g of every
# generator coefficient must pass rho(g) E_Y in Y on its own.  tb4-tight (wave16) tested ~350 single terms at block 4,
# all of depth 1-3 failed, and concluded "multi-term coefficients"; with Theorem T that escape is closed.  Here the
# depth-1 (and optionally depth-2) positions of the ENGINE game are enumerated exhaustively (they are frames of the dual
# polar graph of Sp(8,2): 30 at distance 1, 280 at distance 2 from the start frame) and each representative is tested
# with tb4-tight's own filter9.e_condition (exact mod p = 33554393; a failure is a certain witness when dem = 86,592,
# and if dem > 86,592 every tight schedule is excluded outright anyway).
#
# Verdict depends only on the position P g (P = monomial Cliffords keep N), and P g <-> frame g^-1 L0, L0 = {(0|c)}.
# Labels (a|c) in F_2^4 x F_2^4.  Letter actions on labels (conjugation U X^a Z^c U^dag ~ X^a' Z^c'):
#   dir z = alpha I + beta X^z :  (a|c) -> (a + (c.z) z | c)          (involution mod 2)
#   perm x -> A x + t          :  (a|c) -> (A a | A^-T c)
#   phase i^q(x)               :  (a|c) -> (a | c + S a),  S_ii = q(e_i)-q(0), S_ij = (q(e_i+e_j)-q(e_i)-q(e_j)+q(0))/2  (mod 2)
import sys, os, time, itertools
T0 = time.time()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "make-E"))
import numpy as np
import filter9 as F

DEPTH_MAX = int(os.environ.get("DEPTH_MAX", "1"))
NE = int(os.environ.get("NE", "2"))

def par(x): return bin(x).count("1") & 1

def mat_of_perm(sig):
    t = sig[0]
    cols = [sig[1 << i] ^ t for i in range(4)]          # A e_i
    def A(a):
        r = 0
        for i in range(4):
            if a >> i & 1: r ^= cols[i]
        return r
    # check affinity
    for x in range(16):
        assert sig[x] == A(x) ^ t, "not affine"
    # A^-T: (A^-T c) . (A a) = c . a  -> build table by solving
    AinvT = {}
    for c in range(16):
        for d in range(16):
            if all(par(d & A(1 << i)) == (c >> i & 1) for i in range(4)):
                AinvT[c] = d; break
    return A, AinvT

def S_of_phase(q):
    S = [[0] * 4 for _ in range(4)]
    for i in range(4):
        S[i][i] = (q[1 << i] - q[0]) % 2
        for j in range(i + 1, 4):
            v = (q[(1 << i) | (1 << j)] - q[1 << i] - q[1 << j] + q[0]) % 4
            assert v % 2 == 0, "phase not Clifford-diagonal"
            S[i][j] = S[j][i] = (v // 2) % 2
    return S

def act(letter, v):
    a, c = v
    kind = letter[0]
    if kind in ('dir', 'dirinv'):
        z = letter[1]
        return (a ^ (z if par(c & z) else 0), c)
    if kind == 'perm':
        A, AinvT = mat_of_perm(letter[1])
        return (A(a), AinvT[c])
    if kind == 'phase':
        S = S_of_phase(letter[1])
        Sa = 0
        for i in range(4):
            if sum(S[i][j] * (a >> j & 1) for j in range(4)) % 2: Sa |= 1 << i
        return (a, c ^ Sa)
    raise ValueError(kind)

def act_inv(letter, v):
    kind = letter[0]
    if kind in ('dir', 'dirinv', 'phase'):
        return act(letter, v)                             # all involutions on labels mod 2
    if kind == 'perm':
        sig = letter[1]; inv = [0] * 16
        for x in range(16): inv[sig[x]] = x
        return act(('perm', inv), v)

def lab(v): return v[0] * 16 + v[1]
def unlab(i): return (i >> 4, i & 15)
_tab = {}
def inv_table(letter):
    key = (letter[0], tuple(letter[1]) if isinstance(letter[1], list) else letter[1])
    if key not in _tab:
        _tab[key] = [lab(act_inv(letter, unlab(i))) for i in range(256)]
    return _tab[key]

L0 = frozenset(lab((0, c)) for c in range(16))
def frame(word):
    """frame of g = word[0] ... word[-1] (word[-1] acts first): g^-1 L0, applying w[0]^-1 first"""
    fr = L0
    for let in word:
        t = inv_table(let)
        fr = frozenset(t[i] for i in fr)
    return fr
def depth(fr): return 4 - (len(fr & L0).bit_length() - 1)

# free generators (right multiplication g -> g b): S_i, CZ_ij, coordinate swaps, CNOTs
def q_S(i): return [(x >> i) & 1 for x in range(16)]
def q_CZ(i, j): return [2 * ((x >> i) & 1) * ((x >> j) & 1) for x in range(16)]
def p_swap(i, j):
    out = []
    for x in range(16):
        bi, bj = (x >> i) & 1, (x >> j) & 1
        y = x & ~((1 << i) | (1 << j)) | (bj << i) | (bi << j)
        out.append(y)
    return out
def p_cnot(i, j):   # x_j ^= x_i
    return [x ^ ((((x >> i) & 1)) << j) for x in range(16)]
FREE = [('phase', q_S(i)) for i in range(4)] + [('phase', q_CZ(i, j)) for i in range(4) for j in range(i + 1, 4)] \
     + [('perm', p_swap(i, j)) for i in range(4) for j in range(i + 1, 4)] \
     + [('perm', p_cnot(i, j)) for i in range(4) for j in range(4) if i != j]

def enumerate_depth(dmax):
    """each depth d is one P-orbit (Bruhat cell of Sp(8,2) w.r.t. the Siegel parabolic): find one base word at depth d
    (previous base + dir z), then BFS its P-orbit under right multiplication by free letters (short words)."""
    from collections import deque
    reps = {L0: []}
    base_w, base_fr = [], L0
    for d in range(1, dmax + 1):
        found = None
        for z in range(1, 16):
            t = inv_table(('dir', z))
            fr0 = frozenset(t[i] for i in base_fr)
            if depth(fr0) == d: found = (base_w + [('dir', z)], fr0); break
        assert found, f"no base at depth {d}"
        new = {found[1]: found[0]}
        dq = deque([found])
        while dq:
            u, fr = dq.popleft()
            for b in FREE:
                tb = inv_table(b)
                f2 = frozenset(tb[i] for i in fr)
                if f2 not in new:
                    new[f2] = u + [b]; dq.append((u + [b], f2))
        assert all(depth(f) == d for f in new)
        reps.update(new)
        base_w, base_fr = found
        print(f"depth {d}: {len(new)} frames, longest word {max(len(w) for w in new.values())} letters  ({time.time()-T0:.0f}s)", flush=True)
    return reps

if __name__ == '__main__':
    reps = enumerate_depth(DEPTH_MAX)
    counts = {}
    for fr in reps: counts[depth(fr)] = counts.get(depth(fr), 0) + 1
    print("frames by depth:", dict(sorted(counts.items())), "(dual polar Sp(8,2): 1, 30, 280, 960, 1024)", flush=True)
    assert frame(F.K_WORD) != L0 and depth(frame(F.K_WORD)) == 4
    rng = np.random.default_rng(20261016)
    E0 = F.load_E0(1)
    es = [F.random_E(E0, rng) for _ in range(NE)]
    print(f"{NE} random vectors of E built ({time.time()-T0:.0f}s)", flush=True)
    # controls
    ctrl = [("identity (must pass)", [], True), ("K = d1 d2 d4 d8 (must pass)", F.K_WORD, True),
            ("K then free letters (must pass)", [('perm', p_cnot(0, 1)), ('phase', q_Z := [2 * (x & 1) for x in range(16)])] + F.K_WORD, True)]
    for name, w, want in ctrl:
        got = F.e_condition([(1, w)], es)
        print(f"control {name}: {'PASS' if got else 'FAIL'}  {'ok' if got == want else 'UNEXPECTED'}", flush=True)
    npass = {}
    SL = os.environ.get("SLICE", "0/1"); si, sn = map(int, SL.split("/"))
    onlyd = int(os.environ.get("ONLY_DEPTH", "0"))
    items = [(fr, w) for fr, w in sorted(reps.items(), key=lambda t: (depth(t[0]), sorted(t[0]))) if depth(fr) and (not onlyd or depth(fr) == onlyd)]
    items = items[si::sn]
    print(f"testing {len(items)} positions (slice {SL}, only depth {onlyd or 'all'})", flush=True)
    for fr, w in items:
        d = depth(fr)
        ok = F.e_condition([(1, w)], es)
        npass.setdefault(d, [0, 0]); npass[d][1] += 1
        if ok:
            npass[d][0] += 1
            print(f"PASS depth {d}: word {w}", flush=True)
    for d in sorted(npass):
        print(f"depth {d}: {npass[d][0]} of {npass[d][1]} engine positions pass the f = 9 termwise filter", flush=True)
    print(f"done ({time.time()-T0:.0f}s)")
