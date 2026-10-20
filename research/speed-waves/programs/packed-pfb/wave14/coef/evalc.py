# evaluate an effective coefficient file under a scheme: legality and the 2^21-point metrics
import os, sys, json, numpy as np
P, sc, s, f = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
os.environ['PROTO'] = P; sys.argv = sys.argv[:1]; import gq
sh, gid = gq.setup(sc); c = np.load(f); raw = gq.raw_of(c, sh)
print(json.dumps(dict(file=os.path.basename(f), proto=P, scheme=sc, s=s, legal=gq.legal(raw, gid, s), symmetric=bool(np.all(c == c[::-1])),
      raw_peak=[int(np.abs(raw[gid == j]).max()) for j in range(gid.max() + 1)], **{k: round(v, 3) for k, v in gq.metrics(c).items()})))
