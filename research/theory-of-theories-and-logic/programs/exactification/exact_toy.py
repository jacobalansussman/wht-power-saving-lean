# meta-border, 14 Oct 2026.  Smallest test of the exactification dichotomy (EXACTIFY.md).
# A "unit" Phi(eps; X) is used recursively: U_L = D exactly (the base case), U_j = Phi(eps_j; U_(j+1)).
# Data are truncated Laurent series in ONE eps (matrix coefficients), so we see exactly which coefficient is wrong.
#   formal unit   Phi_f(eps; X) = eps^-1 [ X (I + eps Y + eps^2 N) - X ] Y^-1  = X + eps X N Y^-1
#                 (its eps^-1 terms cancel as a formal identity in X: valid for every X)
#   relational    Phi_r(eps; X) = X + eps N + eps^-1 (X X - S) Z
#                 (its eps^-1 term vanishes only because D D = S, a relation of the paid letter)
# Both have valuation drop a = 1 and are border units at depth 1: Phi(eps; D) = D + O(eps).
# Prediction (EXACTIFY.md s3): with one shared eps, the formal unit's eps^0 coefficient is D at every depth, and
# precision 1 + j*a at depth j suffices; the relational unit's eps^0 coefficient is wrong from depth 2 on,
# whatever the precision; with separated scales eps_j = eps^(c_j) it needs c_(j+1) >= (a+1) c_j.
import numpy as np

ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
S = np.array([[0, 1], [1, 0]], complex)
D = ca * np.eye(2) + cb * S
assert np.allclose(D @ D, S)
rng = np.random.default_rng(3)
Y = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
N = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
Z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
Yi = np.linalg.inv(Y)

# Laurent series: dict {power: 2x2 matrix}, truncated to powers < prec
def smul(A, B, prec):
    out = {}
    for p, a in A.items():
        for q, b in B.items():
            if p + q < prec: out[p + q] = out.get(p + q, 0) + a @ b
    return out
def sadd(*As):
    out = {}
    for A in As:
        for p, a in A.items(): out[p] = out.get(p, 0) + a
    return out
def sc(A, c, shift=0):
    return {p + shift: c * a for p, a in A.items()}
def const(M): return {0: M}

def phi_f(X, prec):   # eps^-1 [ X (I + eps Y + eps^2 N) - X ] Y^-1 ; X known mod eps^(prec+1)
    inner = {0: np.eye(2, dtype=complex), 1: Y, 2: N}
    t = sadd(smul(X, inner, prec + 1), sc(X, -1))
    t = sc(t, 1, shift=-1)
    return {p: a for p, a in smul(t, const(Yi), prec).items() if p < prec}
def phi_r(X, prec):   # X + eps N + eps^-1 (X X - S) Z
    XX = smul(X, X, prec + 1)
    t = sc(smul(sadd(XX, const(-S)), const(Z), prec + 1), 1, shift=-1)
    return {p: a for p, a in sadd(X, {1: N}, t).items() if p < prec}

def run(phi, L, prec_of_depth):
    U = const(D)
    for j in range(L - 1, -1, -1):
        U = phi(U, prec_of_depth(j))
    return U

print("depth L | formal: |coef0 - D| with precision 1+j (enough) | with precision 1 at every depth | relational: |coef0 - D| with precision 1+j, and with precision 1+3j")
for L in range(1, 7):
    ef = np.abs(run(phi_f, L, lambda j: 1 + j).get(0, 0) - D).max()
    ef1 = np.abs(run(phi_f, L, lambda j: 1).get(0, 0) - D).max()
    er = np.abs(run(phi_r, L, lambda j: 1 + j).get(0, 0) - D).max()
    er3 = np.abs(run(phi_r, L, lambda j: 1 + 3 * j).get(0, 0) - D).max()
    print(f"  {L}     | {ef:.2e} | {ef1:.2e} | {er:.2e} | {er3:.2e}")

# separated scales for the relational unit: evaluate numerically at eps_j = e^(c_j), floats
def num_phi_r(X, e): return X + e * N + (X @ X - S) @ Z / e
def num_phi_f(X, e): return (X @ (np.eye(2) + e * Y + e * e * N) - X) @ Yi / e
print("numerical, floats: |U_0 - D| for depth L, shared eps = 1e-3 (formal | relational) and relational with separated scales")
for L in range(1, 5):
    e = 1e-3
    Uf = D.copy(); Ur = D.copy(); Us = D.copy()
    for j in range(L - 1, -1, -1):
        Uf = num_phi_f(Uf, e); Ur = num_phi_r(Ur, e)
        Us = num_phi_r(Us, 1e-2 ** (2 ** j))     # needs c_(j+1) = 2 c_j for a = 1
    print(f"  L={L}: formal {np.abs(Uf - D).max():.2e} | relational shared {np.abs(Ur - D).max():.2e} | relational separated (eps_j = 1e-2^(2^j): needs 2^L digits) {np.abs(Us - D).max():.2e}")
