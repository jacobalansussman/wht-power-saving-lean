#!/bin/bash
# k-build: derive the challenge file for exponent 1 - NUM/10^10 with name suffix SUF from OpenAI's original statement
# (c-base's recipe, BASE-README section CHECKS 3).  usage: derive_challenges.sh <OpenAI ref .lean> <NUM> <SUF> > out.lean
perl -Mutf8 -CSD -pe "s{1 - 1/\(10\^\(13:ℕ\)\)}{1 - $2/(10^(10:ℕ))}; s/\b(decimalExponent|decimalTime|TimeBounds|DFTGoal|ConvGoal|transform_main|convolution_main)\b/\$1$3/g" "$1"
