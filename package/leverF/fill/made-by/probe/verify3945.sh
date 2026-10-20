#!/bin/sh
# o-model: #3945's own owner/verify.py (plain and -O) on a copy of owner/ whose inputs/old-FILL.json.gz and FILL.json.gz are
# a new design's (#3815-stage fill and its #3945 offset regeneration). Only those two files, PROFILE.json (if it changes),
# their SOURCES.json pins and MANIFEST.json's fill hashes (and MANIFEST's own pin) are replaced. Args: NAME OLDFILL NEWFILL NEWPROFILE
set -e
B=/content/work/reopt-o-model; PY=$B/v3815/bin/python; NAME=$1; OLDF=$2; NEWF=$3; NEWP=$4
D=$B/c3945v-$NAME; rm -rf $D; mkdir -p $D; cp -r $B/hub-3945/research/staged544-terminal-gap-fill/owner $D/owner
cp $OLDF $D/owner/inputs/old-FILL.json.gz; cp $NEWF $D/owner/FILL.json.gz; cp $NEWP $D/owner/PROFILE.json
$PY -I -B - $D/owner <<'PYEOF'
import json, gzip, hashlib, sys
from pathlib import Path
R = Path(sys.argv[1]); sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
m = json.loads((R / 'MANIFEST.json').read_text())
m['fill_sha256'] = sha(R / 'FILL.json.gz'); m['fill_payload_sha256'] = hashlib.sha256(gzip.decompress((R / 'FILL.json.gz').read_bytes())).hexdigest()
(R / 'MANIFEST.json').write_text(json.dumps(m, indent=2, sort_keys=True))
s = json.loads((R / 'SOURCES.json').read_text())
for rel in ('inputs/old-FILL.json.gz', 'FILL.json.gz', 'PROFILE.json', 'MANIFEST.json'):
    s['files'][rel] = sha(R / rel)
(R / 'SOURCES.json').write_text(json.dumps(s, indent=2))
print('repinned', {k: s['files'][k][:12] for k in ('inputs/old-FILL.json.gz', 'FILL.json.gz', 'PROFILE.json', 'MANIFEST.json')})
PYEOF
cd $D/owner
nice -n 19 $PY -B verify.py
nice -n 19 $PY -O -B verify.py
echo DONE_VERIFY3945 $NAME
