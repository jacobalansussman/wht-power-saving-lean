#!/bin/sh
# synth.sh LANEFILE SLICEFILE TOP FAMILY(xc7|xcu) -> Yosys synth_xilinx, SYNTHESIS ONLY (as ap-4's synth_xcu.sh, generic).
export PATH=$OSS_CAD_SUITE/bin:$PATH
R=$(cd $(dirname $0)/../rtl; pwd); O=$(cd $(dirname $0); pwd); T=$3; F=$4
LIB="read_verilog -lib +/xilinx/cells_xtra.v"
yosys -q -l $O/e_${T}_${F}.log -p "$LIB; read_verilog $R/$1 $R/$2; hierarchy -top $T; proc; flatten; opt_clean; write_rtlil $O/${T}_${F}.il" || exit 1
yosys -q -l $O/y_${T}_${F}.log -p "$LIB; read_rtlil $O/${T}_${F}.il; synth_xilinx -family $F -noiopad -noclkbuf -top $T; tee -o $O/stat_${T}_${F}.txt stat"
