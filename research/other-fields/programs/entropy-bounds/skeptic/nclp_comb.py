# run nc-gain's nclp.py on a second known-gain network (combination network C(3,2): coding 9, routing 10)
import sys; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'lp-and-layouts'))
import nclp
comb = [('s', 'a'), ('s', 'b'), ('s', 'c'), ('a', 'tab'), ('b', 'tab'), ('a', 'tac'), ('c', 'tac'), ('b', 'tbc'), ('c', 'tbc')]
net = nclp.Net(['s', 'a', 'b', 'c', 'tab', 'tac', 'tbc'], comb, [('s', ['tab', 'tac', 'tbc'], 1), ('s', ['tab', 'tac', 'tbc'], 1)])
lp = nclp.build_lp(net, verbose=False); r = nclp.solve(lp)
print("nclp.py on C(3,2):", r.fun, "(expect 9; routing 10)")
