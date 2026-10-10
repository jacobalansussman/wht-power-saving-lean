"""hio.py -- the data of the search as JSON: hierarchies of shared sums and pools of sums.  Standard library only.

A HIERARCHY (format "e8-hierarchy/1") is the list of shared sums and which target takes which:
   content[u]   what sum u holds: pairs [source, coefficient]            (u < 120: the single source u)
   frame[u]     the span of the labels of its sources (echelon basis, integers = bit vectors of length 9)
   kids[u]      null for a source, else [a, b, r]: sum u = sum a + r * sum b
   need[S]      what target S takes: pairs [sum, coefficient]
   twin         pairs [u, w]: u = a + r b and w = a - r b come from the same two parts (candidates for an in-place pair)
Coefficients are written as strings ("1", "-1/2").  The ORDER of every list is part of the data: the programs walk
these lists in order, and another order can give another circuit.

A POOL (format "e8-pool/1") is a list of sums, each [mask, neg] as two hexadecimal strings: bit T of mask is set when
source T is in the sum, bit T of neg when its coefficient is -1 (the lowest source always has +1).

Files ending in .gz are gzip-compressed, written with a zero time stamp (equal contents give equal bytes)."""
import gzip, json
from fractions import Fraction as Fr

def _read(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        return json.load(fh)

def _write(text, path):
    data = text.encode()
    if path.endswith(".gz"):
        with open(path, "wb") as fh:
            with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0) as gz:
                gz.write(data)
    else:
        with open(path, "wb") as fh:
            fh.write(data)

def load_h(path):
    o = _read(path)
    if o.get("format") != "e8-hierarchy/1":
        raise ValueError("%s: not a hierarchy file (format %r)" % (path, o.get("format")))
    return dict(content=[{int(k): Fr(v) for k, v in c} for c in o["content"]],
                frame=[tuple(int(x) for x in f) for f in o["frame"]],
                kids=[None if k is None else (int(k[0]), int(k[1]), Fr(k[2])) for k in o["kids"]],
                need=[{int(k): Fr(v) for k, v in n} for n in o["need"]],
                twin={int(a): int(b) for a, b in o["twin"]})

def dump_h(d, path, note=None):
    for k in d["kids"]:
        assert k is None or len(k) == 3, "kids must be null or (a, b, r)"
    parts = ['"format":"e8-hierarchy/1"']
    if note: parts.append('"note":' + json.dumps(note))
    js = lambda x: json.dumps(x, separators=(",", ":"))
    parts.append('"content":' + js([[[k, str(v)] for k, v in c.items()] for c in d["content"]]))
    parts.append('"frame":' + js([list(f) for f in d["frame"]]))
    parts.append('"kids":' + js([None if k is None else [k[0], k[1], str(k[2])] for k in d["kids"]]))
    parts.append('"need":' + js([[[k, str(v)] for k, v in n.items()] for n in d["need"]]))
    parts.append('"twin":' + js([[a, b] for a, b in d["twin"].items()]))
    _write("{" + ",\n".join(parts) + "}\n", path)

def load_pool(path):
    o = _read(path)
    if o.get("format") != "e8-pool/1":
        raise ValueError("%s: not a pool file (format %r)" % (path, o.get("format")))
    return [(int(m, 16), int(n, 16)) for m, n in o["sums"]]

def dump_pool(pool, path, note=None):
    parts = ['"format":"e8-pool/1"']
    if note: parts.append('"note":' + json.dumps(note))
    parts.append('"sums":[\n' + ",\n".join('["%x","%x"]' % (m, n) for m, n in pool) + "]")
    _write("{" + ",\n".join(parts) + "}\n", path)
