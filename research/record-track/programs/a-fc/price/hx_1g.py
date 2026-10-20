# hx-price: Theorem 1G constants (PROOF.md section 3) for the full-cut rules: theta, eta_m, kappa, exit bound; [computed]
import json, math, sys
rows = []
for label, f in (('A-fc', 'out/price-Afc.json'), ('design B', 'out/price-controls.json'), ('design A', 'out/price-controls.json')):
    for d in json.load(open(f)):
        if not d['label'].startswith(label): continue
        fc = d['full_cut']; ch = {int(k): v for k, v in fc['children'].items()}; tg = {int(k): v for k, v in fc['targets'].items()}
        b = float(fc['b_routeA_grid24']); tau = 1 - b; m = max(tg); kmin = min(tg); rmax = max(ch)
        eta = min((1 - (1 - t) ** tau) / t ** tau for t in [1 / m + i * (1 - 1 / m) / 200000 for i in range(200001)])
        out = dict(label=d['label'], tau=tau, m=m, k_min=kmin, r_max=rmax, separated=rmax < kmin, paid_rounds=sum(ch.values()),
                   W=sum(tg.values()), theta=math.log(rmax) / math.log(kmin), theta_cap24=math.log(24) / math.log(kmin), eta_m=eta)
        out['kappa'] = math.ceil(math.log(2 / eta) / math.log(1 / out['theta']))
        out['kappa_cap24'] = math.ceil(math.log(2 / eta) / math.log(1 / out['theta_cap24']))
        out['hub_engine_halving_m2_gt_2rmax2'] = m * m > 2 * rmax * rmax
        out['same_with_kmin_in_place_of_m'] = kmin * kmin > 2 * rmax * rmax
        # moment identity at the root (Phi-neutral unit): sum n r^tau == sum W k^tau (float check)
        out['moment_rel_residual_float'] = (sum(n * r ** tau for r, n in ch.items()) - sum(n * k ** tau for k, n in tg.items())) / sum(n * k ** tau for k, n in tg.items())
        rows.append(out)
print(json.dumps(rows, indent=1))
