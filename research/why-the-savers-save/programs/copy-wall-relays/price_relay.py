"""price_relay.py (d-relay) -- WHAT-IF price of the relay payoff the brief names, for scale only.  The relay
decoder it prices does NOT exist in the frame game (no-gain lemma in nogain.py; rank-8 SAT runs in out/).
Engine: m-anatomy levers (Bc for the bit side c of #589, BC = 288 calls; control c0 1.7334e-3; b0 1.6644e-3)."""
import sys, io, contextlib
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').environ.get('PRICER_DIR', '.'))  # folder holding pricer.py (it ships with the lever package's tools)
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../pricing-levers'))  # levers.py
with contextlib.redirect_stdout(io.StringIO()):
    import levers as L
c0, b0 = L.c0, L.b0
print("control c0 = %.10e (bit side #589)   control b0 = %.10e (complex side #799)" % (float(c0), float(b0)))
for name, d in [("P111 cst 63 -> 62 (one rank unit per call)", {7: -1, 6: 1}),
                ("P111 cst 63 -> 56 (relay rank-8 decoder, spans 7: one rank-7 child removed per call)", {7: -1})]:
    c = L.Bc(L.add(L.HB, d, L.BC))
    print("%-88s c = %.4e (%+.2f %%)  what-if, UNREACHABLE by relays" % (name, float(c), float(c / c0 - 1) * 100))
print("kappa: b binds (b0 %.4e < c0 %.4e); a gain in c alone moves kappa by 0; c already sits %+.2f %% above b."
      % (float(b0), float(c0), float(c0 / b0 - 1) * 100))

# Paid-relay variant (relay through a priced descent of T's register; hub C.7: a rank-r descent + re-climb costs
# 2r data-mass units in every call).  Modelled here as 2 extra rank-1 units per call per rank-1 descent; calibration
# line first (hub: one rank-1 descent per local word = -5.630 % c on #555).
one = L.Bc(L.add(L.HB, {1: 2}, L.BC))
print("calibration: one rank-1 descent per call as {1:+2}: c = %.4e (%+.2f %%)  (hub C.7: -5.630 %% c)" % (float(one), float(one / c0 - 1) * 100))
for nd in (4, 17):
    c = L.Bc(L.add(L.HB, {7: -1, 1: 2 * nd}, L.BC))
    print("P111 rank 8 with one rank-7 child removed and %2d rank-1 descents per call: c = %.4e (%+.2f %%)  what-if (>= 17 needed)" % (nd, float(c), float(c / c0 - 1) * 100))
