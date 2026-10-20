"""ctl.py -- controls for the eps-quotient test: (a) two private chains at depth 2 (4 moves) then targets: must PASS;
(b) Theorem A prefixes u3 = D(r u1 + r' u2): must FAIL every target pair (Theorem A, o3-close)."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapeR import *
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
n = int(sys.argv[2]) if len(sys.argv) > 2 else 5
W = World(seed); G = W.G
def geo_letter(M, want):
    for _ in range(3000):
        b = W.rand_letter()
        if W.st(G.mul(W.D, b, M)) == want:
            return b
    return None
t0 = time.time()
y1, y2 = gen_inp(0), gen_inp(1)
for k in range(n):
    u1 = move(W, G.I, y1); u2 = move(W, G.I, y2)
    r = geo_letter(G.D1, (2, 2)); rp = geo_letter(G.D1, (2, 2))
    # (a) private chains
    a = move(W, r, u1); b = move(W, rp, u2)
    res, n1, n2 = test_prefix(W, [y1, y2, u1, u2, a, b])
    print("ctl-private", k, "targets", n1, n2, "pairs", len(res), "pass", sum(x["pass"] for x in res),
          "d4 hist", sorted(collections.Counter(x["d4"] for x in res).items()), flush=True)
    # (b) Theorem A
    u3 = move(W, G.I, comb((1, lmul(G, r, u1, G.I)), (1, lmul(G, rp, u2, G.I))))
    res, n1, n2 = test_prefix(W, [y1, y2, u1, u2, u3])
    print("thmA", k, "targets", n1, n2, "pairs", len(res), "pass", sum(x["pass"] for x in res),
          "d4", sorted(collections.Counter(x["d4"] for x in res).items()),
          "d5", sorted(collections.Counter(x["d5"] for x in res).items()), flush=True)
print("time", round(time.time() - t0, 1))
