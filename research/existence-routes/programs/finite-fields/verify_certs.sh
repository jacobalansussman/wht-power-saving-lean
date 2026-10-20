#!/bin/bash
# verify_certs.sh DIR : for every DIR/<skel>__<sys>.txt, check the certificate exactly against data/dec/<skel>/<sys>
cd "$(dirname "$0")"
for f in "$1"/*.txt; do
  b=$(basename "$f" .txt); skel=${b%%__*}; sys=${b#*__}
  [ -s "$f" ] || { echo "$skel $sys: EMPTY OUTPUT (timeout?)"; continue; }
  printf "%s %s: " "$skel" "$sys"; python3 lift_cert.py check "${DECROOT:-data/dec}/$skel/$sys" "$f" 2>&1 | tail -1
done
