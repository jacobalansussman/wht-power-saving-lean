#!/usr/bin/env python3
"""m-gen: parse the GENERATED Lean data module back (x-mutable's lean_raw_to_gcert.py, text parse only) and check
that it is the input certificate field by field; then run the exact kernel-term mirror (zmirror.py) on the
PARSED data, so that the predicted `decide +kernel` values are those of the text that Lean will read.
usage: roundtrip.py <S0.lean> <cert.gcert.json.gz> <scratch.json.gz>"""
import sys, os, json, gzip, contextlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm
import lean_raw_to_gcert as lr

def main(s0, cert, tmp):
    with contextlib.redirect_stdout(sys.stderr):
        lr.main(s0, tmp)
    p = json.loads(gzip.open(tmp, "rt").read())
    c = json.loads(gzip.open(cert, "rt").read())
    keys = ["h", "v", "R", "ports", "frames", "start", "final", "ext", "ret", "scat", "A", "B"]
    diff = [k for k in keys if p[k] != c[k]]
    out = dict(parsed_equals_certificate=not diff, differing=diff)
    lc = zm.lean_data(p); lc["ext"] = p["ext"]
    s = zm.suite(lc, (10, 21, 4, 8, 8), 4, (10, 8, 4, 8, 8), 80, True)
    out["kernel_terms_on_parsed_text"] = dict(all_true=s["all_true"], scal=[x["value"] for x in s["scal"]],
                                              ychk=s["ychk"]["value"], echk=[x["value"] for x in s["echk"]])
    print(json.dumps(out))
    assert not diff and s["all_true"]

if __name__ == "__main__":
    main(*sys.argv[1:4])
