"""fs-price: independent audit of the GENERATED Lean data of a hosted instance (parses the Lean text; shares no code with
the generator or with host_price).  Rules (h-L3-design / x-hosting, DaysSky's hub #442): 8 runs x 9 stages, call k =
9 run + stage; host i = (child c * 2 + side) * 112 + t, child c in A2 B2 A3 B3 A4 B4 F1 F2 F3 with ranks 16 8 32 8 48 8
48 32 16, pairs 2 2 3 3 4 4 1 2 3, windows (stages) [0,1] [2] [0..3] [4] [0..5] [6] [3..8] [5..8] [7..8]; stage pairs
1 1 1 2 2 3 3 4 4.  SAME-RUN: the host table is one 9 x R table used by every run, so a host only ever holds residuals of
its own run.  Checks: every hosted (stage, slot) lies in its child's window and not in the host's own pair; one slot per
host per stage (role injective); load of each host <= rank and = rank - rem; a slot is hosted at stage l iff it is not
banked in any call k = l (mod 9); one residual per call per bank; bank load = 73 - pad; then the histogram from scratch
(72 x invocation profile from InstU8XM's invPrice, 8 x rem, pads), deficit 2560, exact price (Shea's certify, Sussman's
plan_block at s = 40, A and A + 1).
usage: python3 -I indep_audit.py DATA.lean SCHED.lean INST.lean"""
import sys, re, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
data, sched, inst = [open(f).read() for f in sys.argv[1:4]]
def L(name, t):
    m = re.search(r'def %s : List (?:\(List Nat\)|Nat) := (\[[^\n]*\])\n' % name, t)
    return eval(m.group(1))
hr = L('hrowsL', data)
br = []
for i in range(9):
    br += L('brows%d' % i, data)
pad = L('hpadL', data); rem = L('hremL', data); res = L('schedRes', sched)
R = len(res); nb = len(pad); NH = 2016
RANK = [16, 8, 32, 8, 48, 8, 48, 32, 16]; CPAIR = [2, 2, 3, 3, 4, 4, 1, 2, 3]
WIN = [range(0, 2), [2], range(0, 4), [4], range(0, 6), [6], range(3, 9), range(5, 9), range(7, 9)]
SPAIR = [1, 1, 1, 2, 2, 3, 3, 4, 4]
assert len(hr) == 9 and len(br) == 72 and all(len(r) == R for r in hr + br) and len(rem) == NH
load = Counter()
for l in range(9):
    hosts = [e for e in hr[l] if e < NH]
    assert len(hosts) == len(set(hosts)), 'two slots on one host in one stage'
    for q, e in enumerate(hr[l]):
        if e < NH:
            c = e // 224
            assert l in WIN[c] and CPAIR[c] != SPAIR[l], ('window', l, q, e)
            load[e] += res[q]
        else:
            assert e == NH + q, 'unhosted marker'
for i in range(NH):
    assert load[i] <= RANK[i // 224] and rem[i] == RANK[i // 224] - load[i], ('host load', i)
bload = Counter()
for k in range(72):
    row = br[k]; l = k % 9
    banks = [b for b in row if b < nb]
    assert len(banks) == len(set(banks)), 'two residuals of one call in one bank'
    for q, b in enumerate(row):
        hosted = hr[l][q] < NH
        if b < nb:
            assert not hosted, 'hosted slot also banked'; bload[b] += res[q]
        else:
            assert hosted and b == nb + q, ('unhosted slot without a bank', k, q)
assert all(bload[b] == 73 - pad[b] and 0 <= pad[b] < 73 for b in range(nb)), 'bank load'
m = re.search(r'def invPrice \(φ : ℕ → ℝ\) : ℝ := ([^\n]*)\n', inst)
prof = Counter({int(r): int(n) for n, r in re.findall(r'(\d+) \* φ (\d+)', m.group(1))})
hist = Counter({r: 72 * n for r, n in prof.items()})
hist.update({r: 8 * n for r, n in Counter(rem).items() if r})
hist.update(p for p in pad if p)
W = 64 * 112 + nb
D = 73 * W - sum(r * n for r, n in hist.items())
assert D == 2560, D
sys.path.insert(0, os.path.join(os.environ.get('BC', '.'), 'm-gen/toolkit/tools/gen'))  # see generators/README.txt
from fractions import Fraction as Q
from exact_moment import certify
import foldlib as fl
rate = certify(dict(hist), 73, W)
u = dict(m=73, W=W, D=D, H=dict(hist))
A = int(Q(rate['root_lower']) * 10 ** 10)
while fl.plan_block(u, A, 10, 40) is None: A -= 1
assert fl.plan_block(u, A + 1, 10, 40) is None
hosted = sum(1 for l in range(9) for e in hr[l] if e < NH)
print('INDEPENDENT AUDIT PASS: R %d, banks %d, W %d, hosted entries per run %d (x8 runs = %d residuals hosted of %d), '
      'invocation profile %s, root %s, kernel %d, DFT %d' % (R, nb, W, hosted, 8 * hosted, 72 * R, dict(sorted(prof.items())),
      rate['root_lower_decimal'][:16], A, A - 1))
