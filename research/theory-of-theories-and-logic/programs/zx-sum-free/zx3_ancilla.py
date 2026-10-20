"""zx3_ancilla.py -- a2-zx (meta2, 14 Oct 2026).  Why the ZX model of our game must forbid states and effects on
data wires.  One-bit teleportation, transversal over f columns: data row d (f wires), ancilla row a (f wires) in
|+>^f (phase-free Z-spiders of arity 1, column-uniform), CZ(d_c, a_c) in every column (free), effect <+|^f on the data
row (column-uniform).  The result on the ancilla row is H^(x)f applied to the data row: a paid move for free.
In the game this needs the non-monomial map prod_c (I + X_c)/2 (a sum of 2^f words), which is not a free word."""
import sys
import numpy as np
import pyzx as zx
from pyzx.utils import VertexType, EdgeType

for f in (1, 2, 3):
    g = zx.Graph()
    ins, outs = [], []
    for c in range(f):
        i = g.add_vertex(VertexType.BOUNDARY, qubit=c, row=0)
        dz = g.add_vertex(VertexType.Z, qubit=c, row=1)          # data wire spider; arity: input, CZ leg (effect <+| = no output)
        az = g.add_vertex(VertexType.Z, qubit=f + c, row=1)      # ancilla |+> (arity-1 Z spider) joined to CZ leg and output
        o = g.add_vertex(VertexType.BOUNDARY, qubit=c, row=2)
        g.add_edge((i, dz)); g.add_edge((dz, az), EdgeType.HADAMARD); g.add_edge((az, o))
        ins.append(i); outs.append(o)
    g.set_inputs(ins); g.set_outputs(outs)
    T = zx.tensor.tensor_to_matrix(zx.tensorfy(g, preserve_scalar=False), f, f)
    H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    Hf = H
    for _ in range(f - 1): Hf = np.kron(Hf, H)
    print("f=%d: teleported diagram equals H^(x)f up to scalar: %s" % (f, zx.compare_tensors(T, Hf, preserve_scalar=False)))
