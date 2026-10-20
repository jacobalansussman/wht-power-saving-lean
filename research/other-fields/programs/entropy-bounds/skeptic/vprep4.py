# vprep4.py -- save my Gamma_4 (vgamma) in nc-gain's dpg indexing for the Colab job
import numpy as np, sys
sys.path.insert(0, '.')
from vcut_tests import Mapped
g = Mapped(4)
np.save('gamma4_edges_dpgidx.npy', np.array(g.edges, np.int32)); np.save('gamma4_D_dpgidx.npy', g.D.astype(np.int8))
print(len(g.edges), g.D.shape)
