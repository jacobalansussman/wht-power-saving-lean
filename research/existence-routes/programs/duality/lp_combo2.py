# ex2-duality lp_combo2.py: same LP as lp_combo.py with a wider-window middle relation R2b = 2[4] <= [1] + [7].
import sys
sys.argv = [sys.argv[0]] + sys.argv[1:]
exec(open("lp_combo.py").read().split("R1 = (")[0])
R1 = ({3: 2}, {1: 5})
R2b = ({4: 2}, {1: 1, 7: 1})
Ns = [int(x) for x in sys.argv[1:]] or [30, 60, 120, 240]
for name, rules in [("R2b alone", [R2b]), ("R1 + R2b", [R1, R2b])]:
    tm = theta_max(rules)
    costs = [lp(N, rules) for N in Ns]
    print(f"{name}: duality exponent max Theta = {tm:.4f}; LP cost {[round(c, 4) for c in costs]}; local slopes "
          f"{[round(log(costs[i + 1] / costs[i]) / log(Ns[i + 1] / Ns[i]), 4) for i in range(len(Ns) - 1)]}", flush=True)
