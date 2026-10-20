import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-relays'))
import patterns
patterns.HUB_P111 = [9, 31, 37, 52, 66, 76, 78, 111, 119]
import relaysat
relaysat.instance = patterns.instance
a = sys.argv[1:]
relaysat.solve(a[0], int(a[1]), int(a[2]), a[3], True, None, None, xonly=True)
