#!/bin/bash
REUSE_CHALLENGE=1 /content/work/v-c-wht/vcwht-fresh.sh /content/work/c-wht/lean/complex588 base WHTC588 WHTC588Neg OAI.PowerSaving.WHT.wht_main_block_C588x OAI.PowerSaving.WHTC588.wht_main_C588x > /content/work/logs/v-c-wht/base.log 2>&1
echo "EXIT base $?" >> /content/work/logs/v-c-wht/base.log
/content/work/v-c-wht/vcwht-fresh.sh /content/work/c-wht-q/lean/complex588 packed WHTC588Q WHTC588QNeg OAI.PowerSaving.WHT.wht_main_block_C588Qx OAI.PowerSaving.WHTC588Q.wht_main_C588Qx > /content/work/logs/v-c-wht/packed.log 2>&1
echo "EXIT packed $?" >> /content/work/logs/v-c-wht/packed.log
