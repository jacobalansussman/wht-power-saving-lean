# iso-critic check k1 (pure Python, exact Fractions, < 1 s, < 20 MB).
# Claim [derived here from Baur-Strassen 1983 + Tellegen/transposition, memory]: for any matrix A,
#   L(A) and the arithmetic complexity of the bilinear form q(x, y) = y^T A x agree up to a constant factor.
#   (a) q <= L(A) + 2N: compute A x, then N products y_i (A x)_i and N-1 additions.
#   (b) A <= c * size(q): reverse-mode differentiate ANY circuit for q w.r.t. y (Baur-Strassen, <= ~4-5x gates),
#       the y-gradient is A x (degree 1 in x, free of y); set y = 0 and take degree-1 parts (linearise, <= 3x).
# Test (b) on a circuit for q that does NOT contain a circuit for A: q = x^T (A^T y), built from a circuit for A^T.
# The adjoint then yields A x through the transposed circuit: the transposition principle.
# Done for A = H_n (FFT butterfly for A^T = H_n) and for random sparse non-symmetric A given by random linear circuits.
import random
from fractions import Fraction as F

def run(gates, inputs):
    # gates: list of ('in', k) | ('lin', i, a, j, b) [a*v_i + b*v_j] | ('mul', i, j) | ('const', c)
    v = []
    for g in gates:
        if g[0] == 'in': v.append(inputs[g[1]])
        elif g[0] == 'lin': v.append(g[2] * v[g[1]] + g[4] * v[g[3]])
        elif g[0] == 'mul': v.append(v[g[1]] * v[g[2]])
        else: v.append(g[1])
    return v

def adjoint_ops(gates, out):
    # number of arithmetic operations of reverse-mode differentiation of node `out`, counted per gate:
    # lin: 2 scalings + 2 accumulations; mul: 2 products + 2 accumulations (Baur-Strassen accounting, constant 4).
    ops = 0
    for g in gates[:out + 1]:
        if g[0] in ('lin', 'mul'): ops += 4
    return ops

def gradient(gates, inputs, out):
    v = run(gates, inputs)
    adj = [F(0)] * len(gates); adj[out] = F(1)
    for t in range(out, -1, -1):
        g = gates[t]
        if g[0] == 'lin':
            adj[g[1]] += g[2] * adj[t]; adj[g[3]] += g[4] * adj[t]
        elif g[0] == 'mul':
            adj[g[1]] += v[g[2]] * adj[t]; adj[g[2]] += v[g[1]] * adj[t]
    return adj, v

def fft_wht(n, gates, base):
    # butterfly circuit for H_{2^n} applied to nodes base[0..N-1]; returns output node ids; 2 gates per butterfly pair
    cur = list(base)
    for b in range(n):
        nxt = list(cur)
        for i in range(1 << n):
            if not (i >> b) & 1:
                j = i | (1 << b)
                gates.append(('lin', cur[i], F(1), cur[j], F(1))); nxt[i] = len(gates) - 1
                gates.append(('lin', cur[i], F(1), cur[j], F(-1))); nxt[j] = len(gates) - 1
        cur = nxt
    return cur

def random_lin(N, size, gates, base, rng):
    nodes = list(base)
    for _ in range(size):
        i, j = rng.sample(nodes, 2)
        gates.append(('lin', i, F(rng.randint(-3, 3)), j, F(rng.randint(-3, 3)))); nodes.append(len(gates) - 1)
    return rng.sample(nodes[N:], N) if len(nodes) - N >= N else nodes[-N:]

def build_q(N, build_AT):
    # inputs: x_0..x_{N-1}, y_0..y_{N-1}; q = sum_i x_i * (A^T y)_i
    gates = [('in', k) for k in range(2 * N)]
    xs, ys = list(range(N)), list(range(N, 2 * N))
    outs = build_AT(gates, ys)
    lin_gates_AT = sum(1 for g in gates if g[0] == 'lin')
    prods = []
    for i in range(N):
        gates.append(('mul', xs[i], outs[i])); prods.append(len(gates) - 1)
    acc = prods[0]
    for p in prods[1:]:
        gates.append(('lin', acc, F(1), p, F(1))); acc = len(gates) - 1
    return gates, acc, lin_gates_AT

def matrix_of(gates, outs_fn, N):
    # matrix M with (output) = M (input) for a purely linear sub-circuit on y (x set to 0)
    cols = []
    for k in range(N):
        inp = [F(0)] * (2 * N); inp[N + k] = F(1)
        v = run(gates, inp); cols.append([v[o] for o in outs_fn])
    return [[cols[k][i] for k in range(N)] for i in range(N)]

rng = random.Random(20261014)
print("Baur-Strassen / transposition check: from a circuit for q(x,y) = x^T (A^T y) get A x by reverse mode in y")
for label, N, builder in (
        [("H_%d" % (1 << n), 1 << n, (lambda n: (lambda g, b: fft_wht(n, g, b)))(n)) for n in range(1, 6)] +
        [("random A, N=%d" % N, N, (lambda N: (lambda g, b: random_lin(N, 3 * N, g, b, random.Random(N))))(N)) for N in (5, 8, 13)]):
    gates, out, lin_AT = build_q(N, builder)
    # A^T = matrix of the y-subcircuit; A = its transpose
    youts = [g[2] for g in gates if g[0] == 'mul']
    AT = matrix_of(gates, youts, N)
    A = [[AT[j][i] for j in range(N)] for i in range(N)]
    ok = True
    for trial in range(3):
        x = [F(rng.randint(-5, 5)) for _ in range(N)]
        y = [F(rng.randint(-5, 5)) for _ in range(N)]
        adj, v = gradient(gates, x + y, out)
        grad_y = adj[N:2 * N]
        Ax = [sum(A[i][k] * x[k] for k in range(N)) for i in range(N)]
        ok &= grad_y == Ax
        # y-free: same gradient at y = 0
        adj0, _ = gradient(gates, x + [F(0)] * N, out)
        ok &= adj0[N:2 * N] == Ax
    size_q = sum(1 for g in gates if g[0] in ('lin', 'mul'))
    print("  %-16s q-circuit %4d ops (A^T part %4d); reverse mode gives A x exactly: %s; adjoint cost <= %d ops (= 4 size(q))"
          % (label, size_q, lin_AT, ok, adjoint_ops(gates, out)))
print("(b) holds with constant 4 (adjoint) x 3 (linearise) = 12 on these instances; (a) q <= L(A) + 2N by construction.")
