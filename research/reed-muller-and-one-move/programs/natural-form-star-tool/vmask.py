import struct, sys
raw = open(sys.argv[1], "rb").read(); nf = struct.unpack_from("<i", raw, 0)[0]; o = 4 + 12 * nf
nv = struct.unpack_from("<i", raw, o)[0]; pk = struct.unpack_from("<%dI" % nv, raw, o + 4)
ks = [int(x) for x in sys.argv[2].split(",")]
vm = sum(1 << v for v in range(nv) if all(not (pk[v] >> k & 1) for k in ks))
print(hex(vm))
