"""Price an explicit bank table for Shea's C588 certificate exactly as his pipeline does:
invoice = 72 x Shea's charged local histogram + 112 x bankPrice377 + one block per nonzero pad; W = 64*112 + banks;
exact rational moment root (exact_moment.certify, independent of Lean) and the Lean kernel numerator A on the
1e-10 grid (largest A whose finite-fill s=40 rational certificate passes foldlib.plan_block, i.e. what
generate_global_lean.py emits); DFT numerator = A - 1.  Also validates the table (Shea's checks + controls).
Usage: price_sched.py out.json sched1.schedule.json.gz [sched2 ...]"""
import sys, os, json, gzip
from collections import Counter
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from route_price import IDLE
from exact_moment import certify
import foldlib as fl
from pack_exact import validate

LOCAL = {1: 2596, 2: 492, 3: 280, 4: 232, 5: 56, 6: 56, 8: 120}   # Shea gen-C588 charged local histogram


def invoice(pads):
    hist = Counter({r: 72 * n for r, n in LOCAL.items()})
    for r, n in IDLE.items():
        hist[r] += n * 112
    hist.update(p for p in pads if p)
    W = 64 * 112 + len(pads)
    assert 73 * W - sum(r * n for r, n in hist.items()) == 2560
    return dict(sorted(hist.items())), W


def kernel_A(hist, W, root_lower):
    u = dict(m=73, W=W, D=2560, H=hist)
    A = int(Q(root_lower) * 10 ** 10)
    while fl.plan_block(u, A, 10, 40) is None:
        A -= 1
    return A


def price(sched_path, check=True):
    s = json.loads(gzip.decompress(open(sched_path, 'rb').read()))
    if check:
        validate(s['rows'], s['pads'], s['widths'])
    hist, W = invoice(s['pads'])
    rate = certify(hist, 73, W)
    A = kernel_A(hist, W, rate['root_lower'])
    return dict(schedule=sched_path, banks=len(s['pads']), W=W, pad_hist={str(k): v for k, v in sorted(Counter(s['pads']).items())},
                profile={str(k): v for k, v in hist.items()}, exact_root_lower=rate['root_lower_decimal'][:17],
                exact_grid=int(Q(rate['root_lower']) * 10 ** 10), kernel_A=A, dft=A - 1, rate=rate)


if __name__ == '__main__':
    res = []
    for p in sys.argv[2:]:
        r = price(p)
        print(f"{p}: banks={r['banks']} W={r['W']} pads={r['pad_hist']} exact root={r['exact_root_lower']} "
              f"exact grid={r['exact_grid']} kernel A={r['kernel_A']} DFT={r['dft']} (gap {r['exact_grid']-r['kernel_A']})", flush=True)
        res.append(r)
    json.dump(res, open(sys.argv[1], 'w'), indent=1)
