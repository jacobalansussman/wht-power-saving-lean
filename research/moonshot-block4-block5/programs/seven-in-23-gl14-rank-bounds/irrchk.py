"""irrchk.py -- commutant dimension of sigma(B) (letters) for the address-type representations: irreducible iff 1."""
import sys, os, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4gen import *
L = Lift(seed=21, psi=1, nletters=600)
lets = [b for (b, s) in L.letters]
def comm_dim(sig, ds, k=25):
    rows = []
    for b in lets[:k]:
        S = sig(b); rows.append(np.kron(S, np.eye(ds)) - np.kron(np.eye(ds), S.T))
    M = np.vstack(rows); s = np.linalg.svd(M, compute_uv=False)
    return int((s < 1e-9 * s[0]).sum())
def sig15p(b):
    v = np.ones(16) / 4; U, _, _ = np.linalg.svd(np.eye(16) - np.outer(v, v)); Q = U[:, :15].T
    return Q @ sig_aff16(b) @ Q.T
for name, sig, ds in [('gl14', make_sig('gl14')[0], 14), ('gl15', sig_gl15, 15), ('aff16', sig_aff16, 16), ('15p', sig15p, 15),
                      ('sq16', make_sig('sq16')[0], 16), ('V', sig_nat, 16)]:
    print(name, ds, 'commutant dim', comm_dim(sig, ds), flush=True)
# is sq16 isomorphic to aff16?  intertwiner space dim
rows = []
for b in lets[:25]:
    A = make_sig('sq16')[0](b); B = sig_aff16(b)
    rows.append(np.kron(A, np.eye(16)) - np.kron(np.eye(16), B.T))
s = np.linalg.svd(np.vstack(rows), compute_uv=False); print('Hom(aff16, sq16) dim', int((s < 1e-9 * s[0]).sum()))
# Pauli action: do Paulis act trivially?
from sp8 import sigma as sp
paulis = [b for b in lets if np.allclose(np.abs(b), np.abs(b)) and sp(b) == sp(np.eye(16, dtype=complex))][:5]
print('Pauli-type letters found:', len(paulis))
