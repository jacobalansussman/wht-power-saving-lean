#!/bin/sh
cd "$(dirname "$0")"
WIDE=3 ${PYTHON:-python3} -I lmilp.py ${1:-900} > lmilp_k8_w3.log 2>&1
