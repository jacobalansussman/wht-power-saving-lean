#!/usr/bin/env python3
"""theta.py: the smallest saving that meets OpenAI's challenge.

OpenAI's statement (ComparatorChallenges/UniformFourier.lean) asks for a running time
    n (log n)^theta (log log n)^(4 - theta),   theta = ln(paperMultiplier) / ln(paperM),
    paperMultiplier = paperM - paperDelta / paperW,  paperM = 10^6,  paperW = 2^71,  paperDelta = 6871402692000000.
A saving s (time n (log n)^(1 - s) up to log log factors) meets it as soon as s > 1 - theta.
Pure Python (decimal arithmetic at 60 digits); prints 1 - theta.
"""
from decimal import Decimal, getcontext

getcontext().prec = 60
M = Decimal(10) ** 6
W = Decimal(2) ** 71
delta = Decimal(6871402692000000)
mult = M - delta / W
theta = mult.ln() / M.ln()
print("paperMultiplier =", mult)
print("theta           =", theta)
print("1 - theta       =", 1 - theta)
print("1 - theta (short) = %.4e" % float(1 - theta))
