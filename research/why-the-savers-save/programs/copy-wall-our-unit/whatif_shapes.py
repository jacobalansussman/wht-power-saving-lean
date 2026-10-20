"""whatif_shapes.py (u-ours) -- WHAT-IFS, no design: #799's complex side b with the eight rank-8 copy children of
every call replaced by the copy shapes that our octad/read/rank counting on P112 (toy.py) does not exclude.
Helpers, idle remainder, m and W unchanged (deficit grows by the removed copy mass).  Engine: m-anatomy levers.root
(exact re-implementation of the hub's moment root; control must print 1.6644002789e-3)."""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').environ.get('PRICER_DIR', '.'))  # folder holding pricer.py (it ships with the lever package's tools)
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../pricing-levers'))  # levers.py
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import levers as L
b0 = L.b0
print("control b0 = %.10e (hub #799: 1.6644002789e-3)" % float(b0))
shapes = [("cst 63, one copy 8->7 (hub lever)", {8: -1, 7: 1}),
          ("cst 56, one copy removed (hub lever)", {8: -1}),
          ("cst 63, 9 totals at dim 7 (open by counting)", {8: -8, 7: 9}),
          ("cst 60, 10 totals at dim 6 (open by counting)", {8: -8, 6: 10}),
          ("cst 56, 8 totals at dim 7 (open by counting; excluded for linear totals)", {8: -8, 7: 8}),
          ("cst 54, 9 totals at dim 6 (open by counting)", {8: -8, 6: 9}),
          ("cst 53, octad floor shape 8 x dim 6 + 1 x dim 5", {8: -8, 6: 8, 5: 1}),
          ("cst 0, copies free (star scope dropped)", {8: -8})]
for name, d in shapes:
    b = L.B(L.add(L.H0, d, L.CALLS))
    print("%-80s b = %.4e  (%+.2f %%)  what-if" % (name, float(b), float(b / b0 - 1) * 100))
