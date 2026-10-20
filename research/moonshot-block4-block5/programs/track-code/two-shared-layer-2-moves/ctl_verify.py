"""control for verify.py: the plain 2 x 4 program (8 moves) as a pickle, then verify.py on it"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact import *
from asym import plain
P8 = plain(1)
L = [rec[0][1] for rec in P8.moves]
mv8 = [[(0, L[0])], [(1, L[0])], [(2, L[1])], [(3, L[1])], [(4, L[2])], [(5, L[2])], [(6, L[3])], [(7, L[3])]]
pickle.dump({'moves': mv8}, open('ctl_plain8.pkl', 'wb'))
