"""xc3-cert3: the 405 Fermat planes as lifts of 405 of the 891 alpha-frames (D-coordinates), and their edge excess.
A lift Y_U (3-dim subspace of Q(omega)^6 reducing mod 2 to U) has delta(U,V) = rank[Y_U; Y_V] - 3 >= dist(U,V).
Here: ranks mod p = 2^31 - 1 (omega_p a cube root of 1); mod-p rank <= rank over Q(omega), so this file only
EXPLORES (an excess computed mod p could be too small); the certificate check (cert_check.py) is exact.
"""
import itertools, json, sys
import numpy as np

P = 2147483647
def find_omega():
    for g in range(2, 100):
        w = pow(g, (P - 1) // 3, P)
        if w != 1: return w
OM = find_omega(); assert OM * OM % P != 1 and pow(OM, 3, P) == 1
F4_TO_TEICH = {0: 0, 1: 1, 2: OM, 3: OM * OM % P}   # alpha -> omega

def rank_mod(M):
    M = [list(r) for r in M]
    R = len(M); C = len(M[0]); r = 0
    for c in range(C):
        piv = next((i for i in range(r, R) if M[i][c] % P), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], P - 2, P)
        for i in range(R):
            if i != r and M[i][c] % P:
                f = M[i][c] * inv % P
                M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r])]
        r += 1
        if r == R: break
    return r

def dec(s):
    v = []
    for _ in range(6): v.append(s % 4); s //= 4
    return tuple(reversed(v))
def enc(v):
    s = 0
    for a in v: s = 4 * s + a
    return s

EXP = [1, 2, 3]; LOG = {1: 0, 2: 1, 3: 2}
def fmul(a, b):
    if a == 0 or b == 0: return 0
    return EXP[(LOG[a] + LOG[b]) % 3]
def span(gens):
    S = {(0,) * 6}
    for g in gens:
        S = S | {tuple(x ^ fmul(c, y) for x, y in zip(s, g)) for s in S for c in range(1, 4)}
    return frozenset(enc(v) for v in S)

def matchings(items):
    if not items: yield []; return
    a = items[0]
    for i in range(1, len(items)):
        b = items[i]
        rest = items[1:i] + items[i + 1:]
        for m in matchings(rest): yield [(a, b)] + m

def fermat_planes():
    """returns list of (rows over F4, rows mod p) ; plane: x_a = -zeta x_b for each pair (a, b)"""
    out = []
    for m in matchings(list(range(6))):
        for zs in itertools.product((0, 1, 2), repeat=3):      # zeta = omega^z
            rowsF, rowsP = [], []
            for (a, b), z in zip(m, zs):
                rF = [0] * 6; rP = [0] * 6
                rF[b] = 1; rF[a] = EXP[z]          # reduction of -omega^z is alpha^z
                rP[b] = 1; rP[a] = (-pow(OM, z, P)) % P
                rowsF.append(rF); rowsP.append(rP)
            out.append((m, zs, rowsF, rowsP))
    return out

if __name__ == "__main__":
    data = json.load(open("frames891.json"))
    frames = [frozenset(S) for S in data["frames_D"]]
    fidx = {S: i for i, S in enumerate(frames)}
    Dm = np.load("dist891.npy")
    FP = fermat_planes()
    print("Fermat planes:", len(FP))
    lift = {}
    for m, zs, rF, rP in FP:
        S = span(rF)
        i = fidx.get(S)
        assert i is not None, "Fermat plane not a frame?"
        assert i not in lift, "two Fermat planes reduce to the same frame"
        lift[i] = rP
    print("distinct frames hit:", len(lift))
    fermat = sorted(lift)
    isF = np.zeros(891, bool); isF[fermat] = True
    # edge census
    E = [(i, j) for i in range(891) for j in range(i + 1, 891) if Dm[i, j] == 1]
    print("edges:", len(E))
    ff = sum(1 for i, j in E if isF[i] and isF[j]); fo = sum(1 for i, j in E if isF[i] != isF[j]); oo = len(E) - ff - fo
    print(f"edges Fermat-Fermat {ff}, Fermat-other {fo}, other-other {oo}")
    print("Fermat neighbours of a Fermat frame:", sorted(set(int((Dm[i] == 1)[isF].sum()) for i in fermat)))
    print("Fermat neighbours of an other frame:", sorted(set(int((Dm[i] == 1)[isF].sum()) for i in range(891) if not isF[i])))
    # excess among Fermat pairs at each distance (mod p)
    from collections import Counter
    cnt = {d: Counter() for d in range(4)}
    for a in range(len(fermat)):
        i = fermat[a]
        for b in range(len(fermat)):
            j = fermat[b]
            if j == i: continue
            d = int(Dm[i, j])
            delta = rank_mod(lift[i] + lift[j]) - 3
            cnt[d][delta] += 1
    for d in range(1, 4):
        print(f"Fermat-Fermat ordered pairs at distance {d}: delta histogram {dict(cnt[d])}")
    json.dump({"fermat_idx": fermat, "lift_modp": {str(i): lift[i] for i in fermat}, "p": P, "omega": OM},
              open("fermat405.json", "w"))
