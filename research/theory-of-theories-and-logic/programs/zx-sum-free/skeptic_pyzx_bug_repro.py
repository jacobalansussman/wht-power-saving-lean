"""v-a2-zx-compute: independent reproduction of a2-zx's claimed PyZX 0.10.7 bug (basic_optimization changes the
linear map of a circuit containing XPhase gates).  Replays a2-zx's RNG to rebuild the failing input, then shrinks it
by greedy gate deletion while the bug persists, and prints a minimal failing gate list."""
import sys, random
sys.path.insert(0, ".")
import pyzx as zx
from pyzx.circuit import Circuit
import zx1_sumfree as Z
print("pyzx", getattr(zx, "__version__", "?"), zx.__file__)
rng = random.Random(20261014)
target = (3, 1, 1, 2)
c = None
for m in (2, 3, 4):
    for W in (1, 2, 3):
        for f in (1, 2, 3):
            n = m * W * f
            if n > 36: continue
            for s in range(3):
                A = Z.free_layer(n, 4 * n, rng); B = Z.free_layer(n, 4 * n, rng)
                if (m, W, f, s) == target and c is None:
                    c = Circuit(n); c.add_circuit(B); c.add_circuit(Z.plain(m, W, f)); c.add_circuit(A)
c2 = Circuit.from_graph(zx.simplify.teleport_reduce(c.to_graph()))
print("teleport_reduce+from_graph == input:", zx.compare_tensors(c, c2, preserve_scalar=False))
cb = c2.to_basic_gates()
print("to_basic_gates == c2:", zx.compare_tensors(c2, cb, preserve_scalar=False))
c3 = zx.optimize.basic_optimization(cb.copy())
print("basic_optimization(to_basic_gates(c2)) == c2:", zx.compare_tensors(c2, c3, preserve_scalar=False))
print("XPhase gates in input to basic_optimization:", sum(1 for g in cb.gates if g.name == "XPhase"))

def fails(gates, nq):
    cc = Circuit(nq)
    for g in gates: cc.add_gate(g)
    try:
        out = zx.optimize.basic_optimization(cc.copy())
    except Exception as e:
        return False
    return not zx.compare_tensors(cc, out, preserve_scalar=False)

gates = list(cb.gates)
assert fails(gates, cb.qubits)
changed = True
while changed:
    changed = False
    for i in range(len(gates)):
        trial = gates[:i] + gates[i + 1:]
        if fails(trial, cb.qubits):
            gates = trial; changed = True; break
print("minimal failing input (%d gates on %d qubits):" % (len(gates), cb.qubits))
for g in gates: print("  ", g)
cc = Circuit(cb.qubits)
for g in gates: cc.add_gate(g)
out = zx.optimize.basic_optimization(cc.copy())
print("basic_optimization output:")
for g in out.gates: print("  ", g)
import numpy as np
print("input matrix:\n", np.round(cc.to_matrix(), 3))
print("output matrix:\n", np.round(out.to_matrix(), 3))
