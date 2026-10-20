#!/bin/sh
cd "$(dirname "$0")"
${PYTHON:-python3} -I lmilp.py ${1:-900} > lmilp_k8.log 2>&1
