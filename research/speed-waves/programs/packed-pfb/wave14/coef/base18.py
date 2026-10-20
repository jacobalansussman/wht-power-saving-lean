# Baseline rows: 18-bit (peak 131,071) and 12-bit plain rounding of a prototype, same metrics as gq.py
import os, sys, json, numpy as np
sys.argv = sys.argv[:1]; import gq
for w in (18, 12):
    c = np.round(gq.h / np.abs(gq.h).max() * ((1 << (w - 1)) - 1)).astype(np.int64)
    print(json.dumps(dict(proto=gq.PROTO, case=f'{w}-bit rounded', **gq.metrics(c))))
print(json.dumps(dict(proto=gq.PROTO, case='float', **gq.metrics(gq.h * 1e6))))
