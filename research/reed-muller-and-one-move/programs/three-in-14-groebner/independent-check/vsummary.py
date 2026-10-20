"""vsummary.py: tally every Colab result file of this folder -> out_vsummary.txt."""
import json, glob, collections
out = []
def rows(pat):
    return [json.loads(l) for f in sorted(glob.glob(f"jobs/remote-out/ssh*/{pat}")) for l in open(f)]
def ms(r):
    return [x for v in r["msolve"].values() for x in v]
for name, pat in [("main (20 lifts/set, msolve 2 primes 5 orders/set)", "res_main_[0-9].jsonl"),
                  ("main, every lift (msolve 65521 every order)", "res_mainliftall_*.jsonl"),
                  ("must fail: controls", "res_mainctl_*.jsonl")]:
    R = rows(pat); ns = collections.Counter(str(r["n"]) for r in R)
    m = [x for r in R for x in ms(r)]
    out.append(f"{name}: {len(R)} instances ({sum(c > 1 for c in ns.values())} duplicated), orders {sum(r['orders'] for r in R)}, "
               f"decided {sum(r['decided'] for r in R)}, 1 in ideal over Q {sum(r['unit'] for r in R)}, Singular time-outs {sum(r['timeout'] for r in R)}, "
               f"lifts exact {sum(r['lifts_ok'] for r in R)} of {sum(r['lifts_checked'] for r in R)}, msolve [1] {sum(x is True for x in m)} of {len(m)}")
R = rows("res_mainpass*2_0.jsonl")
for r in R:
    m = ms(r); out.append(f"must pass {r['n']}: orders {r['orders']}, msolve mod 65521 complete non-unit bases {sum(str(x).startswith('BASIS') for x in m)} of {len(m)}")
R = rows("res_mainpass_0.jsonl") + rows("res_mainpassall_0.jsonl")
for r in R:
    out.append(f"must pass {r['n']}: orders {r['orders']}, Singular over Q NOT-UNIT {r['decided'] - r['unit']} of {r['orders']} {sorted(set(r['nonunit'].values()))}")
P = rows("res_plain_0.jsonl") + rows("res_plains_*.jsonl")
out.append(f"must pass plain 3 in 15 (U+, 6 moves), Singular over Q: {len(P)} sets, {sum(r['decided'] for r in P)} orders decided, "
           f"{sum(r['unit'] for r in P)} with 1 in the ideal, time-outs {sum(r['timeout'] for r in P)}; dims {sorted(collections.Counter(v for r in P for v in r['nonunit'].values()).items())}; sets {sorted(r['n'] for r in P)}")
P = rows("res_plainms[01]_*.jsonl")
m = [x for r in P for x in ms(r)]
out.append(f"must pass plain 3 in 15 (U+), msolve mod 65521 (strict wrapper): {len(P)} sets, {len(m)} orders, complete non-unit bases "
           f"{sum(str(x).startswith('BASIS') for x in m)}, [1] {sum(x is True for x in m)}, other {[x for x in m if not (x is True or str(x).startswith('BASIS'))]}; sets {sorted(r['n'] for r in P)}")
open("out_vsummary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
