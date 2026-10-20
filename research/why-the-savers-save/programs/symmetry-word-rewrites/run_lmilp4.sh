#!/bin/sh
cd "$(dirname "$0")"
WIDE=4 ${PYTHON:-python3} -I lmilp.py ${1:-900} > lmilp_k8_w4.log 2>&1
