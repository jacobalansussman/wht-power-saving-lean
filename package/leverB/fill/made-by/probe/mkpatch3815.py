# o-model probe: hall_lp.py of #3815 + a head bonus (our derivative; hub files untouched).
# The bins of class k = 9 r + 4 (children 8/10 of owner r, window (-1, tau5]) that stay empty are literal heads; with the
# half cut each is worth delta(a) = 48^t - 24^t - 73^t + 49^t (t = 1 - a) in moment units, i.e. delta/a in g_cost units.
# Usage: mkpatch3815.py SRC_hall_lp.py DST_hall_lp.py BONUS_SCALE  (1 = the true half-cut value, 0 = control)
import sys
src = open(sys.argv[1]).read(); scale = float(sys.argv[3])
a = "PENALTY = 1e4\n"
b = ("PENALTY = 1e4\n"
     "import math as _m\n"
     "HEAD_SCALE = %r\n"
     "def head_bonus(a):\n"
     "    t = 1 - a\n"
     "    return HEAD_SCALE * (48 ** t - 24 ** t - 73 ** t + 49 ** t) / a\n" % scale)
assert src.count(a) == 1; src = src.replace(a, b)
a = "            costs.append(float(g_cost(left[j], self.a)))\n"
b = ("            costs.append(float(g_cost(left[j], self.a)) - (head_bonus(self.a) if (k % 9 == 4 and int(left[j]) == 48) else 0.0))\n")
assert src.count(a) == 1; src = src.replace(a, b)
# reduced costs in price() must see the same bonus
a = "            rc = g[None, :] - mu[ks, None] + pi[ks] @ P.T.astype(float)\n"
b = ("            rc = g[None, :] - mu[ks, None] + pi[ks] @ P.T.astype(float)\n"
     "            if c == 48:\n"
     "                rc[(ks % 9 == 4)[:, None] & (left == 48)[None, :]] -= head_bonus(self.a)\n")
assert src.count(a) == 1; src = src.replace(a, b)
open(sys.argv[2], 'w').write(src)
