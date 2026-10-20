# Combinatorial facts used by the star lemma, brute force on AG(n,2), n = 3..6:
#  (G1) for every d-flat Psi (d <= n-2) the (d+1)-flats through Psi number 2^(n-d) - 1 and partition AG(n,2) \ Psi;
#  (G2) every set of d+1 points lies in some d-flat.
import itertools
def flats(n, d):
    N = 1 << n; out = set()
    for basis in itertools.combinations(range(1, N), d):
        span = {0}
        for b in basis: span |= {x ^ b for x in span}
        if len(span) != 1 << d: continue
        for t in range(N): out.add(frozenset(x ^ t for x in span))
    return out
for n in range(3, 6):
    N = 1 << n; ok1 = ok2 = True
    fl = {d: flats(n, d) for d in range(0, n)}
    for d in range(0, n - 1):
        for P in fl[d]:
            up = [Phi for Phi in fl[d + 1] if P <= Phi]
            parts = [Phi - P for Phi in up]
            if len(up) != 2 ** (n - d) - 1 or sum(len(q) for q in parts) != N - len(P) or len(set().union(*parts)) != N - len(P):
                ok1 = False
        if d >= 1 and n <= 4:
            for S in itertools.combinations(range(N), d + 1):
                if not any(set(S) <= P for P in fl[d]): ok2 = False; break
    print(f"n={n}: flats per dim {[len(fl[d]) for d in range(n)]}; (G1) {ok1}; (G2) {ok2 if n <= 4 else 'skipped'}")
