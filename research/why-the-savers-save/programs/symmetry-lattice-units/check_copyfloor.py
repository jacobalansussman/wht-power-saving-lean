"""Independent check of the ovoid decoder with the u-floor checker's own Unit methods (copyfloor.py: delivered, fits,
cost; star = rank over Q of the labels read, pairing = rank of the Gram block).  Only the Gram matrix is computed with
numpy (the pure-Python constructor would take hours at v = 2160); everything else is copyfloor's code."""
import sys, os, numpy as np
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker')))
import copyfloor as cf
d = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bitdecoder_terms.npz'))
R, Gm, L = d['R'], d['Gm'], d['L']
U = cf.Unit.__new__(cf.Unit)
U.name, U.field, U.labels, U.p, U.note = 'bw16-bit', 'Q', L.tolist(), 2, ''
U.h = 16; U.v = len(L); U.N = [[int(i == j) for j in range(16)] for i in range(16)]
U.G = (L @ L.T).tolist()
U.knows = [[S != T and U.G[S][T] == 0 for T in range(U.v)] for S in range(U.v)]
U._star, U._pair = {}, {}
terms = [(r.tolist(), g.tolist()) for r, g in zip(R, Gm)]
B = U.delivered(terms)
bad = U.fits(B)
st, pr, rows = U.cost(terms)
print('v', U.v, 'terms', len(terms), 'fits', not bad, 'violations', len(bad), 'star', st, 'pairing', pr, 'x_star', st / U.v)
