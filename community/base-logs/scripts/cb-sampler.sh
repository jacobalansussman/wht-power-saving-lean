#!/bin/bash
# c-base memory sampler: every 5 s, summed RSS (MB) of lean/lake/comparator/lean4export/nanoda processes,
# the largest single one, their count, and free -m "used"/"available". Usage: cb-sampler.sh <log>
LOG=$1
echo "# epoch utc n_proc sum_rss_mb max_rss_mb max_proc mem_used_mb mem_avail_mb" >> $LOG
while true; do
  read n sum mx name < <(ps -eo rss=,comm= | awk '$2 ~ /^(lean|lake|comparator|lean4export|nanoda_bin)$/ {n++; s+=$1; if ($1>m) {m=$1; c=$2}} END {printf "%d %d %d %s\n", n, s/1024, m/1024, (c==""?"-":c)}')
  read used avail < <(free -m | awk '/^Mem:/ {print $3, $7}')
  echo "$(date +%s) $(date -u +%H:%M:%S) $n $sum $mx $name $used $avail" >> $LOG
  sleep 5
done
