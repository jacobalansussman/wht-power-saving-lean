#!/usr/bin/env python3
"""t_strict.py: strict (exact, modulo scalars) coherence of the t_known cases; reuses t_known's word builders by
re-running t_known.main with report() patched to also print the strict class."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ledger, t_known
orig = t_known.report
def rep(name, E, W, R, word, rng, ngauges=4, start=None, deliver=None, out=None):
    g, why = ledger.find_strict_gauge(E, W, R, word, start=start)
    gc, whyc = ledger.find_coherent_gauge(E, W, R, word, start=start)
    res = ledger.run(E, W, R, word, start=start, deliver=deliver)
    print("%-28s m=%d W=%d R=%d M=%2d unit=%s | Sp-class: %-12s | strict: %s" % (name, E.m, W, R, res["M"], res["unit"],
          whyc.split(" (")[0], why))
    if out is not None: out.append(dict(case=name, rows=[dict(raw_ok=True, unit_ok=True, unit=res["unit"])]))
    return []
t_known.report = rep
t_known.main()
