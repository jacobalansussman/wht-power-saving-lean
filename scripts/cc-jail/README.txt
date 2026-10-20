cc-jail (prep, wht45/lever, 19 Oct 2026): the stand-in jail for the OFFICIAL comparator runs, as fs-build / fsp-cc /
fsp-build used it (NOT the official landrun/systemd sandbox): namespaces, read-only mounts except <project>/.lake and a
private /tmp (JAIL_TMP), uid 65534, no capabilities, no_new_privs, socket(AF_UNIX) blocked by nounix (seccomp).
  jail2.sh       = fs-build's jail2.sh, only the nounix path changed to /content/work/cc-jail/nounix
  nounix.c       sha256 f90dbc3996f3c56b... (same source as before); nounix compiled here with gcc -O2 (hash in HASHES)
  fresh-bin/lake sha256 01f56ebe57de64f1... (LEAN_NUM_THREADS=1 wrapper of the v4.34.1 lake)
Recipe: fsp-build's colab/fspb-cmp.sh (kept in my private working notes, not in this repository), with
F=/content/work/<label> and these paths: it chowns YOUR COPY's .lake to 65534 for the run (chown back to root before
any root lake), and needs /content/work/cc/bin (comparator, lean4export, nanoda_bin) from BASE-READY.
