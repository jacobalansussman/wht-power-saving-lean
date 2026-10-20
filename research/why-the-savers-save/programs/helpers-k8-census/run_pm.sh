#!/bin/sh
# k8-pairs: d-pair-market's pm.py (copied unchanged) on the k8 demand (k8_demand.json == a0_demand.json here) and U8 partition.
cd "$(dirname "$0")"
VP=${PYTHON:-python3}
for B in 0 1 2 3 4 5 6 7 8 9 10 11 12 13; do
  env OPTS=P,LL,LF,H FB=8 FACELOW=0 TL=60 DUMP=own/b${B}_sol.json $VP -I pm.py $B > own/b$B.json
done
for B in 2 7; do for FB in 0 4; do
  env OPTS=P,LL,LF,H FB=$FB FACELOW=0 TL=60 $VP -I pm.py $B > own/b${B}_fb$FB.json
done; done
echo done > own/DONE
