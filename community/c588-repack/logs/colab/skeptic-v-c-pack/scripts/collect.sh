#!/bin/sh
# v-c-pack: copy the remaining logs (compare-neg) back from the Colab.
H=~/misc/fft/verify/checks/wht43/colab-run
$H/ssh.sh 'tail -6 /content/work/logs/v-c-pack/vcp2.log; grep -E "^FAIL|RESULT" /content/work/logs/v-c-pack/compare-neg.log'
$H/get.sh '/content/work/logs/v-c-pack/*' $H/logs/v-c-pack/
