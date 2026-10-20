#!/bin/sh
# workdir.sh [DEST]: build a working copy of this section's programs, with each folder's data next to its programs.
#
# The programs read and write their input and output files in the current folder, and a few read files of a
# sibling folder (for example ../natural-form-solver/). So this script makes one folder, DEST (default: work/ in
# this section's folder), with programs/<name>/ and data/<name>/ merged into DEST/<name>/, and unpacks every
# job-outputs.tar.gz into a jobs/ folder beside it and systems.tar.gz into gbin/. Run the commands of each
# folder's README from inside DEST/<name>/. Nothing outside DEST is changed.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
SEC=$(cd "$HERE/.." && pwd)
DEST=${1:-$SEC/work}
mkdir -p "$DEST"
cp -R "$HERE/." "$DEST/"
rm -f "$DEST/workdir.sh"
cp -R "$SEC/data/." "$DEST/"
find "$DEST" -name 'job-outputs.tar.gz' | while read -r t; do
  d=$(dirname "$t"); mkdir -p "$d/jobs"; tar xzf "$t" -C "$d/jobs"; rm -f "$t"
done
find "$DEST" -name 'systems.tar.gz' | while read -r t; do
  d=$(dirname "$t"); tar xzf "$t" -C "$d"; rm -f "$t"
done
find "$DEST" -name '*.sh' -exec chmod +x {} \;
echo "working copy: $DEST"
