"""fs-price: a standalone Lean module that runs XHost.hostSchedK (h-L3-build, soundness HostSched.ofTab_fitsR) on the
generated U8XMHH tables with every list inline (no certificate import), plus a negative twin (one remainder off by one).
usage: python3 -I mk_mac_test.py DATA.lean SCHED.lean OUTDIR"""
import sys, re, os
data, sched, out = open(sys.argv[1]).read(), open(sys.argv[2]).read(), sys.argv[3]
def L(name, t):
    return re.search(r'def %s : List (?:\(List Nat\)|Nat) := (\[[^\n]*\])\n' % name, t).group(1)
hr = L('hrowsL', data); pad = L('hpadL', data); rem = L('hremL', data); res = L('schedRes', sched)
br = ' ++ '.join('brows%d' % i for i in range(9))
nb = len(eval(pad)); R = len(eval(res))
remv = eval(rem); i = next(j for j, x in enumerate(remv) if x > 0); bad = list(remv); bad[i] -= 1
for name, rr, want in (('HsU8XMHH', rem, 'true'), ('HsU8XMHHNeg', str(bad), 'false')):
    t = 'import Work.XHost.HSchedK\n/-! fs-price: hostSchedK on the generated U8XMHH tables (lists inline). -/\n'
    t += 'set_option maxRecDepth 100000\nset_option maxHeartbeats 4000000\nnamespace FsPrice.%s\nopen OAI.PowerSaving.XHost\n' % name
    t += 'noncomputable def hrowsL : List (List Nat) := %s\n' % hr
    for k in range(9):
        t += 'noncomputable def brows%d : List (List Nat) := %s\n' % (k, L('brows%d' % k, data))
    t += 'noncomputable def browsL : List (List Nat) := %s\n' % br
    t += 'noncomputable def hpadL : List Nat := %s\nnoncomputable def hremL : List Nat := %s\nnoncomputable def res : List Nat := %s\n' % (pad, rr, res)
    t += 'theorem t_ok : hostSchedK %d %d hrowsL browsL hpadL hremL res = %s := by decide +kernel\nend FsPrice.%s\n' % (R, nb, want, name)
    open(os.path.join(out, name + '.lean'), 'w').write(t)
print('written', R, nb, 'negative: rem[%d] %d -> %d' % (i, remv[i], bad[i]))
