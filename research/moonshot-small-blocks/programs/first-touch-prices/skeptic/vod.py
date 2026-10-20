# vod.py -- v-w7-mixed3: own computation of the over-demand OD, the Lemma Q first-touch bound and the Lemma Q'' flag
# bound for small programs in the hook (level k = 0), to check w7-mixed3's out_lemmaq_check.txt and the slack identity
# cap*M - sum dem = sum (cap - gain) + OD.
import sys, numpy as np
sys.argv = [sys.argv[0], 'engine', '1', '0']
exec(open(__file__.rsplit('/', 1)[0] + '/vhook.py').read().split('game = sys.argv[1]')[0])
P = len(NI)
def R(g): return rho(g)[np.ix_(NI, NI)], rho(g)[np.ix_(OI, NI)]   # (Y -> Y part, Y -> Y' part) restricted to N
I8 = np.eye(8)
def run(name, W, targets, moves):
    # moves: list of dicts {slot: 8x8 matrix} = slot coefficients c_{t,l} of generator t (already expanded)
    dimY = W * P
    # Kb: v in N^W with pi sum_l rho(c_tl) v_l = 0 for all t
    rows = []; gains = []; cur = np.zeros((0, dimY), complex)
    for mv in moves:
        blk = np.zeros((len(OI), dimY), complex)
        for l, g in mv.items(): blk[:, l * P:(l + 1) * P] = R(g)[1]
        new = np.vstack([cur, blk]); gains.append(np.linalg.matrix_rank(new, tol=1e-7) - np.linalg.matrix_rank(cur, tol=1e-7) if len(cur) else np.linalg.matrix_rank(blk, tol=1e-7)); cur = new
    dimKb = dimY - np.linalg.matrix_rank(cur, tol=1e-7)
    Es = [E[t] for t in targets]; dimE = sum(e.shape[1] for e in Es)
    dem = sum(P - e.shape[1] for e in Es)
    # Kb inside (+)E ?  check: Kb basis projected
    OD = dimE - dimKb
    s = cap * len(moves) - dem
    # Lemma Q'': U_1 = N^W; r_t = rank pi C_t on U_t cap (+)E; U_{t+1} = U_t cap ker C_t (full action, X-valued)
    Ebig = np.zeros((dimY, dimE), complex); c0 = 0
    for w, e in enumerate(Es): Ebig[w * P:(w + 1) * P, c0:c0 + e.shape[1]] = e; c0 += e.shape[1]
    U = np.eye(dimY, dtype=complex); rq = []
    for mv in moves:
        if U.shape[1] == 0: rq.append(0); continue
        full = np.zeros((560, dimY), complex); proj = np.zeros((len(OI), dimY), complex)
        for l, g in mv.items():
            Rg = rho(g); full[:, l * P:(l + 1) * P] = Rg[:, NI]; proj[:, l * P:(l + 1) * P] = Rg[np.ix_(OI, NI)]
        # basis of U cap (+)E
        M = np.hstack([U, -Ebig]); ns = nullsp(M, 1e-8); UE = U @ ns[:U.shape[1]] if ns.shape[1] else np.zeros((dimY, 0))
        rq.append(int(np.linalg.matrix_rank(proj @ UE, tol=1e-7)) if UE.shape[1] else 0)
        FU = full @ U; ns2 = nullsp(FU, 1e-8) if np.abs(FU).max() > 1e-12 else np.eye(U.shape[1]); U = U @ ns2
    print(f'{name}: gains {gains}, demand {dem}, cap*M - dem = {s}, sum(cap-gain) = {sum(cap - g for g in gains)}, OD = {OD}, '
          f'identity holds: {s == sum(cap - g for g in gains) + OD}; Lemma Q\'\' {sum(rq)} {rq}; untouched dim {U.shape[1]}', flush=True)
Dd = {1: D1, 2: D2, 4: D4}
RR = D1 @ D2
# plain 2 R: u1 = D1 x1, u2 = D2 u1, u3 = D1 x2, u4 = D2 u3
run('plain 2R', 2, ['R', 'R'], [{0: D1}, {0: D2 @ D1}, {1: D1}, {1: D2 @ D1}])
# packed: u1 = D1(x1 + x2), u2 = D1 x1, u3 = D2 u2, u4 = D2 (u1 - u2)
run('packed 2R', 2, ['R', 'R'], [{0: D1, 1: D1}, {0: D1}, {0: D2 @ D1}, {1: D2 @ D1}])
# plain K: D1 x, D2 u1, D4 u2
run('plain K', 1, ['K'], [{0: D1}, {0: D2 @ D1}, {0: D4 @ D2 @ D1}])
# a hook-impossible control: 1 K + 1 R in 4 moves (shared): u1 = D1(x1 + x2) ... just check OD > s detection with an R in 1 move
