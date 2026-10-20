"""price.py (d-bitspan) -- WHAT-IFS on PR589's bit side c: the nine rank-7 copy children of every call ({7: 9} of the
per-call copy entry, 288 calls) replaced by other copy shapes.  Engine: m-anatomy levers (Bc, HB, add, BC); the
control must print c0 = 1.7334270e-3.  No construction is priced here: every row is a what-if."""
import sys, io, contextlib; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').environ.get('PRICER_DIR', '.'))  # folder holding pricer.py (it ships with the lever package's tools)
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../pricing-levers'))  # levers.py
with contextlib.redirect_stdout(io.StringIO()):
    import levers as Lv
c0 = Lv.Bc(Lv.HB); b0 = Lv.b0
print("control c0 = %.7e (PR589 bit side, expected 1.7334270e-3); b0 = %.7e (#799 complex side)" % (float(c0), float(b0)))
shapes = [("cst 62: one copy 7 -> 6 (lever check, expect +2.80 %)", {7: -1, 6: 1}),
          ("cst 62: 8 totals at 6 + 2 at 7 (first door, K = 10)", {7: -7, 6: 8}),
          ("cst 62: 9 totals at 6 + 1 at 8 (K = 10)", {7: -9, 6: 9, 8: 1}),
          ("cst 61: 9 totals at 6 + 1 at 7 (K = 10)", {7: -8, 6: 9}),
          ("cst 60: 10 totals at 6 (K = 10)", {7: -9, 6: 10}),
          ("cst 62: 7 at 6 + 4 at 5 (K = 11)", {7: -9, 6: 7, 5: 4}),
          ("cst 56: 8 totals at 7 (clique-count floor shape)", {7: -1}),
          ("cst 0: copies free", {7: -9})]
for name, d in shapes:
    c = Lv.Bc(Lv.add(Lv.HB, d, Lv.BC))
    print("%-62s c = %.4e  (%+.2f %%)  what-if" % (name, float(c), float(c / c0 - 1) * 100))
print("kappa: b0 = %.4e < c0 = %.4e, so b binds; a gain in c moves kappa only together with a complex-side gain above "
      "%+.2f %% on b." % (float(b0), float(c0), float(c0 / b0 - 1) * 100))
print("P120 bit decoder (v = 120, a different unit, not a profile delta): x = cst/v; break-even with P111's 63/111 = %.4f "
      "is cst = %.1f; today's P120 floor 72 gives x = %.4f (estimate only)." % (63 / 111, 120 * 63 / 111, 72 / 120))
