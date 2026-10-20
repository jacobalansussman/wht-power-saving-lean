import nclp
arcs = [('s1', 't1'), ('s1', 'a'), ('s2', 'a'), ('s2', 't2'), ('a', 'b'), ('b', 't1'), ('b', 't2')]
net = nclp.Net(['s1', 's2', 'a', 'b', 't1', 't2'], arcs, [('s1', ['t1', 't2'], 1), ('s2', ['t1', 't2'], 1)])
r = nclp.solve(nclp.build_lp(net))
res, canon, key, active = nclp.solve_cp(net, batch=500)
print("butterfly multicast: full LP", r.fun, " cutting-plane", res.fun)
