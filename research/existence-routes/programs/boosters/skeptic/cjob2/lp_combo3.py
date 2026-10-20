# ex2-duality lp_combo3.py: LP of lp_combo.py for the spreading relation Rs = 2[2] <= [1] + [3] with and without R1.
import sys
exec(open("lp_combo.py").read().split("R1 = (")[0])
R1 = ({3: 2}, {1: 5})
Rs = ({2: 2}, {1: 1, 3: 1})
Ns = [int(x) for x in sys.argv[1:]] or [30, 60, 120, 240]
for name, rules in [("Rs alone", [Rs]), ("R1 + Rs", [R1, Rs])]:
    tm = theta_max(rules)
    costs = [lp(N, rules) for N in Ns]
    print(f"{name}: duality exponent max Theta = {tm:.4f}; LP cost {[round(c, 4) for c in costs]}; local slopes "
          f"{[round(log(costs[i + 1] / costs[i]) / log(Ns[i + 1] / Ns[i]), 4) for i in range(len(Ns) - 1)]}", flush=True)
