"""Lemma O check, independent: |tr((c K c')^* w)|/N for one-paid-letter words w and frames c, c' (dense, f = 2..5)."""
import sys, numpy as np
from mini import Game
frames = ["1", "CS", "A", "P0", "lCS", "chi0", "chi1i", "C01.CS", "P1.chi0"]
def mat(G, fn):
    return fn(np.eye(G.N, dtype=complex))
for f in range(2, 6):
    G = Game(f)
    Km = mat(G, G.K)
    # unitarity / flatness of K
    assert np.allclose(Km.conj().T @ Km, np.eye(G.N)); assert np.allclose(np.abs(Km), G.N ** -0.5)
    tops = []
    for c in frames:
        for c2 in frames:
            tops.append(mat(G, lambda X: G.word(c, G.K(G.word(c2, X)))))
    ws = {}
    for z in (1, 2, 3):
        for u in ["1", "CS", "A", "P0", "chi1"]:
            for v in ["1", "CS", "P1", "chi0i"]:
                ws["%s.D%d.%s" % (u, z, v)] = mat(G, lambda X: G.word(u, G.D(z, G.word(v, X))))
    free = {w: mat(G, lambda X: G.word(w, X)) for w in ["1", "CS", "A", "P0", "chi0.CS", "P0.P1"]}
    m1 = max(abs(np.vdot(T, Wm)) / G.N for T in tops for Wm in ws.values())
    m0 = max(abs(np.vdot(T, Wm)) / G.N for T in tops for Wm in free.values())
    # two-paid words (not covered by Lemma O): K itself, D1 P1 D2 CS
    two = {"D1.D2": mat(G, lambda X: G.D(1, G.D(2, X))),
           "D1.P1.D2.CS": mat(G, lambda X: G.D(1, G.word("P1", G.D(2, G.word("CS", X))))),
           "D1.CS.D2": mat(G, lambda X: G.D(1, G.word("CS", G.D(2, X))))}
    mt = {k: max(abs(np.vdot(T, Wm)) / G.N for T in tops) for k, Wm in two.items()}
    print("f=%d  max one-paid %.4e (2^-f/2 = %.4e)  max free %.4e (2^-f = %.4e)  two-paid %s" %
          (f, m1, 2 ** (-f / 2), m0, 2.0 ** -f, " ".join("%s:%.4f" % kv for kv in mt.items())), flush=True)
