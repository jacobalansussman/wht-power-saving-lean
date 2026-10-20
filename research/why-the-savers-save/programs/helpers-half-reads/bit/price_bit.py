"""d-halfread BIT SIDE pricing (what-if bracket; the hub's hosting pipeline was not rerun): c = certify(288 x local
paid profile + idle remainder REM, 64, 28416, bad envelope) with the hub's bit moment engine (copied as bit_moment.py
from research/e8-bit-carrier-skeleton/.../engines/moment.py), exactly as #589's tune/fastc.py does.  Control: #589's c.
Retiring carrier 470 moves 8 units per call into the idle remainder; where they land depends on the re-run hosting, so
every landing (a fresh rank-8 child, or a bin leftover L in REM growing to L+8) is priced.  Usage: python3 -I price_bit.py"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bit_moment as M
REM = {2: 13920, 4: 14208, 6: 14208, 7: 21312, 42: 288}; T = 288; m = 64; W = 28416
base = {1: 1323, 2: 1342, 3: 79, 4: 44, 5: 102, 6: 11, 7: 19}
new = json.load(open(os.path.join(HERE, 'retire470.json')))['new']; new = {int(k): v for k, v in new.items()}
def glob(loc, rem):
    P = {}
    for r, n in loc.items(): P[r] = P.get(r, 0) + T*n
    for r, n in rem.items(): P[r] = P.get(r, 0) + n
    P = {r: n for r, n in P.items() if n}
    assert m*W - sum(r*n for r, n in P.items()) == 10272; return P
c0 = M.certify(glob(base, REM), m, W, True)['lower']
print('control c', c0, ' #589 tune/EXPECTED: 0.001733427039589844234749  equal:', c0 == '0.001733427039589844234749')
res = {'control': c0}
for L in [0, 2, 4, 6, 7, 42]:
    rem = dict(REM)
    if L: rem[L] -= T
    rem[L + 8] = rem.get(L + 8, 0) + T
    c = M.certify(glob(new, rem), m, W, True)['lower']; res['leftover %d -> %d' % (L, L + 8)] = c
    print('freed 8 units/call land on a bin leftover %2d -> %2d: c = %s (%+.4f %% of c)' % (L, L + 8, c, (float(c)/float(c0) - 1)*100))
json.dump(res, open(os.path.join(HERE, 'price_bit.json'), 'w'), indent=1)
