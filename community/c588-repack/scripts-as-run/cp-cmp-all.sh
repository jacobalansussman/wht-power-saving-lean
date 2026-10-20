#!/bin/bash
# c-pack: official comparator (stand-in sandbox jail2.sh) on the built overlay tree: positive (Shea's config unchanged)
# then negative (challenge one unit stronger). .lake handed to uid 65534 first, as c-base did for t-c588.
P=/content/work/cpack/lean/complex588
s=$(date +%s); chown -R 65534:65534 $P/.lake; echo "chown wall=$(( $(date +%s)-s ))s"
/content/work/cp-cmp.sh P1-fourier-shea-config $P comparator/fourier.json; echo "P1 rc=$?"
/content/work/cp-cmp.sh N1-negative $P /content/work/cc/configs/fourier-neg.json; echo "N1 rc=$? (must be nonzero)"
