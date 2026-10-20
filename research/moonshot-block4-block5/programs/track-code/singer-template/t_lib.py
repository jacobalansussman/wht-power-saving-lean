"""t_lib.py: calibration of s4lib (Singer lift order, frame keys, plain words found by the joint solver, negative control)."""
import time, sys
from s4lib import *
t0 = time.time()
S = S4(4, seed=2); G = S.G
print("eta lift order (projective):", S.eta_order_proj, " eta symplectic P =", S.P, f"({time.time()-t0:.1f}s)")
for _ in range(100):
    g = G.mul(S.D1, S.rand_free(), S.D1, S.rand_free())
    b = S.rand_free()
    assert fast_frame_key(G.mul(b, g)) == fast_frame_key(g)
print("frame key invariance ok")
e = S.eta; ei = np.conj(e).T
print("eta D1 eta^-1 in B' D1 (eta fixes N<x>):", fast_frame_key(G.mul(S.D1, ei)) == fast_frame_key(G.mul(e, S.D1, ei)),
      "| K^-1 eta^-1 K monomial (eta in H_max):", G.is_monomial(G.mul(np.conj(G.K).T, ei, G.K)))
def coord(perm):
    sig = [sum((((x >> j) & 1) << perm[j]) for j in range(4)) for x in range(16)]
    return G.mono(sig, [0] * 16)
def chain(l, perms):
    g = S.inp(l); out = []
    for p in perms:
        g = S.move([(1, coord(p), g)]); out.append(g)
    return out
perms = [(0, 1, 2, 3), (1, 0, 2, 3), (2, 1, 0, 3), (3, 1, 2, 0)]
gens = chain(0, perms[:3]) + chain(1, perms[:3])
print("gen depths/dK:", [(sorted(set(S.depth(gid) for (l, gid) in g)), sorted(set(S.dK(gid) for (l, gid) in g))) for g in gens])
pool = E_pool(S, gens)
print('E pool size', len(pool))
Es = [pool[0], pool[1 if len(pool) > 1 else 0]]
for name, alpha in (("identity", [[{0: 1}, {}], [{}, {0: 1}]]),
                    ("halving", [[{0: Fr(1, 2)}, {0: Fr(1, 2)}], [{0: Fr(1, 2)}, {0: Fr(-1, 2)}]])):
    t1 = time.time()
    r, info = S.joint(gens, 2, Es, alpha)
    print("plain control,", name, "->", r, (info if r != 'yes' else info['st']), f"{time.time()-t1:.2f}s", flush=True)
t1 = time.time()
r, info = S.joint(gens[:3] + gens[3:5], 2, Es, [[{0: 1}, {}], [{}, {0: 1}]])
print("negative control (array 2 one move short):", r, info if r != 'yes' else info['st'], f"{time.time()-t1:.2f}s")
