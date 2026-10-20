# vprep.py -- save my Gamma_m (vgamma) in nc-gain's dpg indexing for Colab jobs
import numpy as np, sys
sys.path.insert(0, '.')
from vcut_tests import Mapped
for m in map(int, sys.argv[1:]):
    g = Mapped(m)
    np.save(f'gamma{m}_edges_dpgidx.npy', np.array(g.edges, np.int32)); np.save(f'gamma{m}_D_dpgidx.npy', g.D.astype(np.int8))
    print(m, len(g.edges), g.D.shape)
