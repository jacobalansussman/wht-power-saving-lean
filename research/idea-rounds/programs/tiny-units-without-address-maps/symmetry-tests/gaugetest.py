#!/usr/bin/python3
"""gaugetest.py (y-theory): exact test of the gauge of S3b on the straight-line form of real schedules.
(1) x_t -> u x_t with u = 2 shift(5);  (2) x_t -> Z x_t with Z = phase(p,0)^2, p.z_t = 1, where dir Z dir^-1 is another free map;
(3) x_t -> x_t + shift(3) x_j for an earlier move j with the same direction.  Then replay by regmachine."""
src = open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "symtest.py")).read().rsplit("\nmain()", 1)[0]
exec(src)
def prod(letters):
    a = np.eye(8, dtype=complex)
    for lt in letters:
        a = rm.letter_matrix(M, (lt[0], 0) + tuple(lt[1:])) @ a
    return a
def lmul(mat, letters, combo, real, c=1):
    out = {}
    for (sym, k), v in combo.items():
        k2 = rm.key(mat @ rm.MATS[k]); real.setdefault(k2, real[k] + letters)
        out[(sym, k2)] = out.get((sym, k2), 0) + v * c
    return out
def rsub(xs, fin, t, real, pieces):
    """replace every reference kappa.y_t by sum over pieces (sym2, mat, letters, c) of c kappa.mat.sym2"""
    def fix(combo):
        out = {}
        for (sym, k), v in combo.items():
            if sym == ("y", t):
                for sym2, mat, letters, c in pieces:
                    k2 = rm.key(rm.MATS[k] @ mat); real.setdefault(k2, letters + real[k])
                    out[(sym2, k2)] = out.get((sym2, k2), 0) + v * c
            else:
                out[(sym, k)] = out.get((sym, k), 0) + v
        return {a: b for a, b in out.items() if b}
    return [fix(x) if i > t else x for i, x in enumerate(xs)], [fix(x) for x in fin]
def replay(dirs, xs, fin, real):
    w2, n2 = rebuild(dirs, xs, fin, real)
    f = rm.formal_run(M, n2, LIVE, w2)
    num = rm.live_kernel_numeric(M, n2, LIVE, w2, fs=(1, 2, 3), trials=1)
    assert all(v == f for v in num.values())
    return f, len(moves(w2))
allok = True
for name, w, n, want in [("plain9", plain9(), 3, True), ("wasteful11", wasteful(), 4, True), ("invalid_mix", invalid_mix(), 4, False), ("invalid8", invalid8(), 3, False)]:
    dirs, xs, fin, real = extract(w, n)
    s = len(dirs)
    # (1) unit 2.shift(5) on a move in the middle
    t = s // 2; b = rm.shift(M, 5)
    xs1 = list(xs); xs1[t] = lmul(b, [("shift", 5)], xs[t], real, c=2)
    d = rm.dirm(M, dirs[t]); bp = d @ b @ np.linalg.inv(d); assert np.array_equal(bp, b)
    xs1, fin1 = rsub(xs1, fin, t, real, [(("y", t), np.linalg.inv(bp), [("shift", 5)], Fr(1, 2))])
    r1 = replay(dirs, xs1, fin1, real)
    # (2) Z = phase(p,0)^2 with p.z_t = 1 on move t = 1: dir Z dir^-1 = i Z shift(z), a different free map
    t = 1; z = dirs[t]; p = [q for q in range(1, 8) if rm.dot(q, z) == 1][0]
    Zm = rm.phase_mat(M, p, 0) @ rm.phase_mat(M, p, 0); d = rm.dirm(M, z); bp = d @ Zm @ np.linalg.inv(d)
    assert not np.array_equal(bp, Zm) and rm.key(bp) is not None
    inv_letters = [("phase", p, 0), ("phase", p, 1)] * 3 + [("phase", p, 0)] * 2 + [("shift", z)]
    assert np.array_equal(prod(inv_letters), np.linalg.inv(bp)), "letters do not give (dir Z dir^-1)^-1"
    xs2 = list(xs); xs2[t] = lmul(Zm, [("phase", p, 0)] * 2, xs[t], real)
    xs2, fin2 = rsub(xs2, fin, t, real, [(("y", t), np.linalg.inv(bp), inv_letters, Fr(1))])
    r2 = replay(dirs, xs2, fin2, real)
    # (3) triangular: x_t += shift(3) x_j, j < t the first earlier move with the same direction
    pair = [(j, t) for t in range(s) for j in range(t) if dirs[j] == dirs[t]][0]
    j, t = pair; tau = rm.shift(M, 3)
    xs3 = list(xs); add = lmul(tau, [("shift", 3)], xs[j], real)
    x = dict(xs[t])
    for kk, v in add.items(): x[kk] = x.get(kk, 0) + v
    xs3[t] = {a: c for a, c in x.items() if c}
    xs3, fin3 = rsub(xs3, fin, t, real, [(("y", t), np.eye(8, dtype=complex), [], Fr(1)), (("y", j), tau, [("shift", 3)], Fr(-1))])
    r3 = replay(dirs, xs3, fin3, real)
    ok = r1[0] == want and r2[0] == want and r3[0] == want and r1[1] == r2[1] == r3[1] == s
    allok &= ok
    print("%s  S3b gauge %s: unit 2.shift -> %s; Z with dir Z dir^-1 = i Z shift (phase vector %d on direction %d) -> %s; x_%d += shift(3) x_%d -> %s; expected %s; moves %d"
          % ("PASS" if ok else "FAIL", name, r1[0], p, z, r2[0], t, j, r3[0], want, s))
print("ALL PASS" if allok else "SOME FAIL")
