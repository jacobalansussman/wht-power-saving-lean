#!/usr/bin/env python3
"""Compare fresh outputs of tools/reproduce-search.sh with the recorded results.

usage: python3 -I -B tools/compare_results.py FRESH_DIR RECORDED_DIR
Both folders have the subfolders heads/ and g0/. A file that parses as JSON is compared as JSON, with every string
that contains a '/' (a file path) replaced by its last component, because the runs used other folders. Other files
are compared byte for byte. Timing files (*.time, *.err) are skipped. Prints one line per file and, at the end,
"ALL EQUAL" or the number of differences (exit code 1).
"""
import json, os, sys

fresh, rec = sys.argv[1], sys.argv[2]


def norm(x):
    if isinstance(x, dict):
        return {k: norm(v) for k, v in x.items()}
    if isinstance(x, list):
        return [norm(v) for v in x]
    if isinstance(x, str) and '/' in x:
        return x.rstrip('/').rsplit('/', 1)[-1]
    return x


def load(p):
    b = open(p, 'rb').read()
    try:
        return 'json', norm(json.loads(b))
    except ValueError:
        return 'bytes', b


bad = 0
for sub in ('heads', 'g0'):
    for name in sorted(os.listdir(os.path.join(rec, sub))):
        if name.endswith(('.time', '.err')):
            continue
        r = os.path.join(rec, sub, name)
        f = os.path.join(fresh, sub, name)
        if not os.path.exists(f):
            print(f'  not re-run   {sub}/{name}')
            continue
        kr, vr = load(r)
        kf, vf = load(f)
        same = kr == kf and vr == vf
        print(f"  {'equal  ' if same else 'DIFFERS'}      {sub}/{name}")
        bad += not same
print('ALL EQUAL' if not bad else f'{bad} DIFFERENCES')
sys.exit(1 if bad else 0)
