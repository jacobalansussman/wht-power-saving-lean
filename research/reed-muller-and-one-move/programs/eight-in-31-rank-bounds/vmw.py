"""vmw.py tables ks t: hex mask of parity classes v with w(v) = #{l : <class v, T_k_l> = 1} <= t."""
import struct, sys
raw = open(sys.argv[1], "rb").read(); nf = struct.unpack_from("<i", raw, 0)[0]; o = 4 + 12 * nf
nv = struct.unpack_from("<i", raw, o)[0]; pk = struct.unpack_from("<%dI" % nv, raw, o + 4)
ks = [int(x) for x in sys.argv[2].split(",")]; t = int(sys.argv[3])
w = [sum(pk[v] >> k & 1 for k in ks) for v in range(nv)]
if len(sys.argv) > 4: print(sorted(w))
print(hex(sum(1 << v for v in range(nv) if w[v] <= t)))
