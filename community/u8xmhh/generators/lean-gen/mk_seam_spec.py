"""fs-price: seam spec for mk_fourier_hx.py (the hosted instance NEW of certificate OLD), numbers from the generator's JSON.
usage: python3 -I mk_seam_spec.py OLD NEW SUMMARY.json OLD_KERNEL OLD_W OUT.json"""
import sys, json
OLD, NEW, summ, oldA, oldW, out = sys.argv[1], sys.argv[2], json.load(open(sys.argv[3])), sys.argv[4], sys.argv[5], sys.argv[6]
A, W = summ['kernel_numerator'], summ['W']
spec = {"keep": {"GXD.%s." % NEW: "GXD.%s." % OLD, "blocks%sx" % NEW: "blocks%sx" % OLD, "fold_%sx" % NEW: "fold_%sx" % OLD},
        "num": {oldA: str(A), str(int(oldA) - 1): str(A - 1), oldW: str(W)},
        "import": "Work.XHost.%s" % NEW,
        "call": ("  XHost.engine_of_ginv_host GXD.{O}.ginv GXD.{O}.h9 GXD.{N}.hsched GXD.{N}.hfits\n"
                 "    73 {W} rfl GXD.{N}.hwidth (by norm_num) BlockAccounting.blocks{N}x\n"
                 "    (fun φ => by rw [GXD.{N}.hprice φ]; exact BridgeRate.{N}.unit_cost_x φ)\n"
                 "    40 (1 - {A}/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_{N}x cl m k v\n").format(O=OLD, N=NEW, W=W, A=A)}
json.dump(spec, open(out, 'w'))
print(json.dumps(spec)[:300])
