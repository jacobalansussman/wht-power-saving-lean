import nclp
arcs = [('s1', 't1'), ('s1', 'a'), ('s2', 'a'), ('s2', 't2'), ('a', 'b'), ('b', 't1'), ('b', 't2')]
net = nclp.Net(['s1', 's2', 'a', 'b', 't1', 't2'], arcs, [('s1', ['t1', 't2'], 1), ('s2', ['t1', 't2'], 1)])
lp = nclp.build_closed_lp(net)
print("butterfly multicast: simplex", nclp.solve(lp, verbose=False).fun, " pdlp (primal, rigorous LB):", nclp.solve_pdlp(lp, 2)[:2])
net = nclp.Net([0, 1, 2, 3], [(0,1),(1,0),(1,2),(2,1),(2,3),(3,2)], [(0, [3], 1)])
lp = nclp.build_closed_lp(net)
print("path:", nclp.solve_pdlp(lp, 1)[:2])
