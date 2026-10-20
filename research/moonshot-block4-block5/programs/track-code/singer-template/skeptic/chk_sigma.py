"""chk_sigma.py: the order-3 Levi element (S4_ORDER=3) used for order-3 filters: order, commutation with D1 and K."""
import os, numpy as np
os.environ['S4_ORDER'] = '3'
from s4lib import S4
S = S4(4, seed=1); G = S.G
c = S.eta
print('P =', S.P, 'matrix order', S.eta_order_proj, 'len(etap)', len(S.etap), 'monomial', G.is_monomial(c))
print('commutes with D1:', np.allclose(G.mul(c, S.D1), G.mul(S.D1, c)))
print('K c K^-1 monomial (c in H_max):', G.is_monomial(G.mul(G.K, c, np.linalg.inv(G.K))))
print('commutes with K:', np.allclose(G.mul(c, G.K), G.mul(G.K, c)))
