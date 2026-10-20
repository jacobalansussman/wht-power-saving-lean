#!/bin/bash
# copy the built base tree, give its .lake to user 65534, run the official comparator on Shea's unchanged config
set -x
s=$(date +%s)
cp -a /content/work/shea /content/work/cc/t-c588 || exit 91
echo "copy wall=$(( $(date +%s)-s ))s"
chown -R 65534:65534 /content/work/cc/t-c588/lean/complex588/.lake
/content/work/cb-cmp.sh A1-fourier-shea-config /content/work/cc/t-c588/lean/complex588 comparator/fourier.json
