#!/bin/bash
# v-c-pack: run a command in /content/work/cpack/lean/complex588 with the whole /content/work/cpack tree bind-mounted read-only (private mount namespace).
TREE=/content/work/cpack; P=$TREE/lean/complex588
exec unshare --mount --propagation private -- /bin/bash -c 'mount --bind '"$TREE $TREE"' && mount -o remount,bind,ro '"$TREE"' && cd '"$P"' && exec "$@"' _ "$@"
