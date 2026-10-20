#!/bin/bash
# c-base: build the official comparator's tools at the pins of bigospeedrun/verifier profiles/tools.json
# (= Shea's receipt/verification.json): comparator d03acab, lean4export 076e8e5, nanoda 3a24072, landrun 811cfff.
# Toolchain file of comparator and lean4export set to v4.34.1 (as VERIFY.md and the verifier's setup() do).
set -u
. /content/work/env.sh
CC=/content/work/cc
mkdir -p $CC/bin $CC/dl $CC/home
cd $CC
echo "== START $(date -u) =="
gclone() { # gclone <url> <sha> <dir>
  git init -q "$3" && git -C "$3" -c core.hooksPath=/dev/null fetch -q --depth=1 "$1" "$2" && git -C "$3" -c core.hooksPath=/dev/null checkout -q --detach FETCH_HEAD
  echo "$3: $(git -C "$3" rev-parse HEAD) (wanted $2)"
}
gclone https://github.com/leanprover/comparator.git d03acab154d269c06e60e4de7e4cc85deebff94b $CC/comparator
gclone https://github.com/leanprover/lean4export.git 076e8e57707e813375e8f9da8bf989799ace9680 $CC/lean4export
gclone https://github.com/ammkrn/nanoda_lib.git 3a2407216ee84a75f9e1aead6803d0578be06ae7 $CC/nanoda
gclone https://github.com/Zouuup/landrun.git 811cfff51ceaf3d9843708aa6d22e9b84ccac8b4 $CC/landrun
echo "== lean4export =="
cd $CC/lean4export; echo "leanprover/lean4:v4.34.1" > lean-toolchain; git diff --stat
/usr/bin/time -v -o $CC/dl/l4e.time lake build 2>&1 | tail -5; grep -E "Maximum res|Elapsed" $CC/dl/l4e.time
echo "== comparator =="
cd $CC/comparator; echo "leanprover/lean4:v4.34.1" > lean-toolchain; git diff --stat
grep -n -A3 lean4export lake-manifest.json | head -8
/usr/bin/time -v -o $CC/dl/cmp.time lake build 2>&1 | tail -5; grep -E "Maximum res|Elapsed" $CC/dl/cmp.time
git -C .lake/packages/lean4export rev-parse HEAD 2>/dev/null
cp -v $CC/lean4export/.lake/build/bin/lean4export $CC/comparator/.lake/build/bin/comparator $CC/bin/
echo "== rust (rustup into /content/work/rust) =="
export RUSTUP_HOME=/content/work/rust/rustup CARGO_HOME=/content/work/rust/cargo
curl -sSf https://sh.rustup.rs -o $CC/dl/rustup-init.sh; sha256sum $CC/dl/rustup-init.sh
sh $CC/dl/rustup-init.sh -y --no-modify-path --profile minimal --default-toolchain stable 2>&1 | tail -3
export PATH=$CARGO_HOME/bin:$PATH
cargo --version; rustc --version
echo "== nanoda =="
cd $CC/nanoda
/usr/bin/time -v -o $CC/dl/nanoda.time cargo build --release --locked -j 4 2>&1 | tail -5; grep -E "Maximum res|Elapsed|Exit" $CC/dl/nanoda.time
cp -v $CC/nanoda/target/release/nanoda_bin $CC/bin/
echo "== go =="
cd $CC/dl
curl -sSfL "https://go.dev/dl/?mode=json" -o go-versions.json
GOV=$(python3 -I -c "import json;d=json.load(open('go-versions.json'));print([r for r in d if r['stable']][0]['version'])")
GOSHA=$(python3 -I -c "import json;d=json.load(open('go-versions.json'));r=[r for r in d if r['stable']][0];print([f for f in r['files'] if f['os']=='linux' and f['arch']=='amd64' and f['kind']=='archive'][0]['sha256'])")
echo "go version chosen: $GOV sha256 $GOSHA"
curl -sSfL "https://go.dev/dl/$GOV.linux-amd64.tar.gz" -o go.tgz
echo "$GOSHA  go.tgz" | sha256sum -c -
mkdir -p $CC/goroot && tar -C $CC/goroot -xzf go.tgz
export GOROOT=$CC/goroot/go GOPATH=$CC/gopath GOCACHE=$CC/gocache GOTOOLCHAIN=local HOME=$CC/home
export PATH=$GOROOT/bin:$PATH
go version
echo "== landrun =="
cd $CC/landrun
/usr/bin/time -v -o $CC/dl/landrun.time go build -o $CC/bin/landrun ./cmd/landrun 2>&1 | tail -5; grep -E "Exit" $CC/dl/landrun.time
echo "== nounix =="
gcc -O2 -Wall -o $CC/bin/nounix $CC/nounix.c && echo nounix built
ls -la $CC/bin; sha256sum $CC/bin/*
$CC/bin/landrun --version
echo "== END $(date -u) =="
