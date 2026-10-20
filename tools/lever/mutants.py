#!/usr/bin/env python3
"""g-heads: mutants that tools/lower_b.py must reject (each one either raises in the lowering / checks, or flips a
reported equality to false). usage: mutants.py CANONICAL-WORD.json.gz BORROW-WORD.json.gz OUT.json"""
import sys, json, copy, importlib.util, os
from collections import Counter
spec = importlib.util.spec_from_file_location('lb', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lower_b.py'))
lb = importlib.util.module_from_spec(spec); spec.loader.exec_module(lb)
c0 = lb.read(sys.argv[1]); word = lb.read(sys.argv[2])
hw = lb.hub_widths(word, c0['helper_origin'])


def verdict(c):
    try:
        cert, _ = lb.lower(c)
        d = lb.uncopy(cert)
        chk = lb.check(d, c0['local_profile'])
        reg = lb.total_matrix(d)
        ec = lb.exit_conditions(d, reg, False)
        dc = lb.dirt_check(c, d)
        res = lb.residuals(d)
    except AssertionError as e:
        return dict(rejected=True, by='assertion: ' + str(e)[:160])
    fails = [k for k, ok in (('hist == hub profile', chk['hist_equals_hub_native_local_profile']),
                             ('R', ec['R_rows_offdiag_zero']), ('C', ec['C_cols_slot_offdiag_zero']),
                             ('hid', ec['hid_y_from_x_identity']), ('hy', ec['hy_y_columns_unit']),
                             ('dirt', dc['y_rows_exactly_y_plus_x']), ('widths', res == hw)) if not ok]
    return dict(rejected=bool(fails), by=fails)


out = {'control_unmutated': verdict(c0)}
assert not out['control_unmutated']['rejected'], out
v = c0['v']
# 1. a B gate's frame replaced by a frame of smaller dimension (non-monotone climb)
m = copy.deepcopy(c0)
i = next(i for i, e in enumerate(m['B']) if e[0] == 'add' and len(m['frames'][e[-1]]) >= 3)
m['B'][i][-1] = next(f for f, fr in enumerate(m['frames']) if len(fr) == 1)
out['B_gate_frame_lowered'] = verdict(m)
# 2. a warp operator's G changed (two columns swapped)
m = copy.deepcopy(c0)
w = m['warps'][0]; g = w['G_columns']; g[3], g[4] = g[4], g[3]
out['warp_G_columns_swapped'] = verdict(m)
# 3. one initial correction dropped (the twisted dirt is then not cancelled)
m = copy.deepcopy(c0)
k = next(i for i, e in enumerate(m['precondition']) if e[0] == 'tadd')
del m['precondition'][k]
out['initial_correction_dropped'] = verdict(m)
# 4. one compensation coefficient doubled
m = copy.deepcopy(c0)
m['compensation'][0][3] *= 2
out['compensation_coefficient_doubled'] = verdict(m)
# 5. one exit helper ends at the full frame instead of its hyperplane (width one too high)
m = copy.deepcopy(c0)
r = m['zero_buffers'][0]
m['final'][r] = 1
out['exit_helper_widened_to_full'] = verdict(m)
# 6. one cut coefficient changed (scatter wrong)
m = copy.deepcopy(c0)
m['cut'][0][3] += 1
out['cut_coefficient_changed'] = verdict(m)
# 7. one lifted-target add moved to a non-lifted y
m = copy.deepcopy(c0)
lifted = {int(t) + v for t in m['target_levels']}
i = next(i for i, e in enumerate(m['A']) if e[0] == 'add' and e[1] in lifted)
m['A'][i][1] = next(r for r in range(v, 2 * v) if r not in lifted)
out['lifted_add_on_unlifted_target'] = verdict(m)
out['all_mutants_rejected'] = all(x['rejected'] for k, x in out.items() if k != 'control_unmutated')
json.dump(out, open(sys.argv[3], 'w'), indent=1)
print(json.dumps(out, indent=1))
