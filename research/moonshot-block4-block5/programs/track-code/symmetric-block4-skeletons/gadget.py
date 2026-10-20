"""gadget.py -- mw-sym4: exhaustive test of the SYMMETRIC TWO-LAYER GADGET (theory in STATUS section E).
H cyclic of order W acting regularly on W arrays through c in H_max.  Layer 0: one tau-fixed shared move
  v = delta * sum_k chi^k tau^k(b0 x_0)            (delta a tau-fixed direction, b0 a progress letter, chi^W = 1)
Layer 1: one free orbit  u^(j) = tau^j(D1 a v)     (a a progress letter for v's x_0 term).
The gadget works at the 2-step frames iff (ii) for every j the x_j-term of u^(0) is free-equivalent to the x_j-term of
u^(j)'s 'own' progressed term (all cross terms coincide), and (iii) the resulting W x W twisted circulant over C[B] is
invertible.  Then the C[B]-span of the layer-1 orbit holds, for every array, its input advanced by TWO steps, at cost
1 + W (instead of 2W): with W more ... private moves per array this gives (m-2)W + W + 1 = (m-1)W + 1 moves.
This script counts, over all (delta, b0, a) and characters chi, how many satisfy (i) progress on x_0 for both layers
and (ii); for those it reports the coincidence pattern and the order of the circulant's generator eta."""
import sys, time
from symsolve import *

def order_of(G, M, cap=64):
    P = M.copy(); I = np.eye(G.n)
    for k in range(1, cap + 1):
        if np.allclose(P, I):
            return k
        P = G.mul(P, M)
    return None

def main():
    m, W, Hn = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    samples = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
    t0 = time.time()
    R = Run(m, W, Hn, 1); G = R.G; rng = random.Random(5)
    # generator h = first nontrivial element whose pi is the cycle l -> l+1 (BFS order: powers of the generator)
    Hc = R.Hc
    fd = fixed_dirs(R, rng, samples)
    tot = prog = coinc = 0
    found = []
    for d in fd:
        b0s = G.progress_letters(G.ID, rng, want=24, Dmat=d, steps=1)
        for b0 in b0s:
            g0 = G.intern(G.mul(d, b0))                     # v's x_0 term
            # v as a dict: sum_k tau^k(b0 x0), coefficient 1 (characters handled below symbolically)
            base = {(0, G.intern(b0)): 1}
            terms = []                                      # (k, l, gid) term of tau^k(b0 x0) after delta
            for k, h in enumerate(Hc):
                img = R.tau(h, base) if k else base
                (l, gb), _ = next(iter(img.items()))
                terms.append((k, l, G.intern(G.mul(d, G.MATS[gb]))))
            avs = G.progress_letters(g0, rng, want=24, Dmat=G.D1, steps=1)
            for a in avs:
                tot += 1
                # u0 = D1 a v : its terms (k, l, gid)
                u0 = [(k, l, G.intern(G.mul(G.D1, a, G.MATS[g]))) for (k, l, g) in terms]
                own = [t for t in u0 if t[1] == 0][0]       # progressed x_0 term
                if G.dist(own[2]) != m - 2:
                    continue
                prog += 1
                # images u^(j) = tau^j(u0): own term of u^(j) is tau^j(own) on input j
                ok = True; etas = []
                for (k, l, g) in u0:
                    if l == 0:
                        continue
                    # find the orbit element h_j with pi(0) = l
                    hj = next(h for h in Hc if h[1][0] == l)
                    lj, gj = R.tau_term_cached(hj, 0, own[2])
                    b = G.mul(G.MATS[g], np.conj(G.MATS[gj]).T)
                    if not G.is_monomial(b):
                        ok = False; break
                    etas.append((l, k, b))
                if not ok:
                    continue
                coinc += 1
                if len(found) < 20:
                    found.append((d, b0, a, etas))
    print(f"m={m} W={W} H={Hn}: fixed dirs {len(fd)}; (delta,b0,a) tested {tot}; two-step progress on x_0 {prog}; "
          f"ALL cross terms coincide (ii) {coinc} ({time.time()-t0:.0f}s)", flush=True)
    for (d, b0, a, etas) in found[:5]:
        print("  coincidence: cross terms (input l from tau^k-component) = beta * own_l ; orders of beta:",
              [(l, k, order_of(G, b)) for (l, k, b) in etas])

if __name__ == "__main__":
    main()
