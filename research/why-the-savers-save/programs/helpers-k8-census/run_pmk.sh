#!/bin/sh
# k8-pairs: K-aware per-block MILP (pmk.py) on the six leaver-cell blocks of k8.
cd "$(dirname "$0")"
VP=${PYTHON:-python3}
env OPTS=P,LL,LF,H FB=8 FACELOW=0 TL=60 $VP -I pmk.py 1 > kaware/ctrl_b1_fb8.json
for B in 0 2 5 7 8 10; do for FB in 0 8; do
  env OPTS=P,LL,LF,H FB=$FB FACELOW=0 TL=90 DROP=cellneeds_k8.json NOLOW=nolow_k8.json DUMP=kaware/b${B}_fb${FB}_sol.json $VP -I pmk.py $B > kaware/b${B}_fb$FB.json
done; done
echo done > kaware/DONE
