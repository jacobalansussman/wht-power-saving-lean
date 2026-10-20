# p-alt: compare the warp sets and initial corrections of BORROW-WORD files (data only).
# Key per warp without its B position: (lender, G columns, frame columns ...). Prints a digest per word, and a digest
# of the (lender, G) pairs alone (lenderG), which is the same on every borrowed word.
import sys, json, gzip, hashlib
def load(p):
    b = open(p, 'rb').read(); return json.loads(gzip.decompress(b) if b[:2] == b'\x1f\x8b' else b)
def dig(x): return hashlib.sha256(json.dumps(x, sort_keys=True).encode()).hexdigest()[:12]
for p in sys.argv[1:]:
    d = load(p)
    w = d['warps']; ic = d['initial_corrections']
    keys = sorted(json.dumps(x[1:]) for x in w)          # drop the B position
    pos = [x[0] for x in w]
    lenders = sorted({x[1] for x in w})
    nontriv = sum(1 for x in w if x[2] != [1 << i for i in range(len(x[2]))])
    lg = sorted(json.dumps(x[1:3]) for x in w)          # (lender, G) only: frames and positions dropped
    print(f"{p.split('/dl/')[-1]:90s} warps={len(w)} nonidentityG={nontriv} lenders={len(lenders)} lenderG={dig(lg)} warpset={dig(keys)} Bpos[{min(pos)}..{max(pos)}] ic={len(ic)} icset={dig(sorted(json.dumps(c, sort_keys=True) for c in ic))} removed={dig(sorted(d['removed']))}")
