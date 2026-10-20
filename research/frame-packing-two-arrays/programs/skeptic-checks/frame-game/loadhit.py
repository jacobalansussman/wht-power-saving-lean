"""loadhit.py: load fp-build's hit pickle with a restricted unpickler (numpy reconstruction globals only)."""
import pickle, sys
import numpy as np
ALLOWED = {("numpy._core.multiarray", "_reconstruct"), ("numpy.core.multiarray", "_reconstruct"),
           ("numpy", "ndarray"), ("numpy", "dtype"), ("numpy._core.multiarray", "scalar"), ("numpy.core.multiarray", "scalar")}
class R(pickle.Unpickler):
    def find_class(self, mod, name):
        if (mod, name) not in ALLOWED: raise pickle.UnpicklingError(f"blocked {mod}.{name}")
        if mod.startswith("numpy._core") and not hasattr(np, "_core"): mod = mod.replace("_core", "core")
        return super().find_class(mod, name)
def load(path):
    with open(path, "rb") as fh: return R(fh).load()
if __name__ == "__main__":
    d = load(sys.argv[1])
    for k, v in d.items():
        if k != "inj": print(k, v)
    ident = tuple(1 << i for i in range(2 * d["m"]))
    c = d["c"]
    for i, (rg, sl, src, mu) in enumerate(d["inj"]):
        print(i, rg, sl, src, "identity chart" if tuple(mu) == ident else mu, np.round(c[i], 4))
