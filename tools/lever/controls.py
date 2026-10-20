# lever-reprice: reproduce the hub's complex roots b for the newest top words, digit for digit.
# Two independent formulas: (1) m-anatomy's exact re-implementation of the hub moment engine (price_b, 1e-24 grid,
# certified by rational enclosures); (2) rl-plan's mixed-target Decimal root (rescue.py's formula, generalised to any W).
# Run: python3 -I -B tools/controls.py
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pricer import price_b, dec
from mixed import root
WORDS = {
    # PR: (profile, W, hub's stated complex root / kappa headline, where stated)
    '1839 T74': ({1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328, 48: 22748}, 66304,
                 '0.0016797549561255445 (rl-plan control)'),
    '3148 T740': ({1: 16376940, 2: 4821840, 3: 1567431, 4: 679320, 5: 147075, 8: 53280, 48: 282569}, 663040,
                  'kappa 0.001729777494678597011886 - 1e-46 - 1e-70 (PR body); profile from evidence/FINITE.json on main'),
    '3605 T740': ({1: 16417896, 2: 4835160, 3: 1540174, 4: 679321, 5: 143992, 6: 3, 8: 53280, 48: 283185}, 663040,
                  'root [0.001730006549538841847005, ...006] (PR body); histogram published in #3659'),
    '3642 T530': ({1: 11747920, 2: 3453480, 3: 1112042, 4: 486540, 5: 103690, 8: 38160, 48: 202828}, 474880,
                  'kappa 0.001730394391813963415120 - 1e-46 - 1e-70 (PR body)'),
    '3653 T488': ({1: 10817496, 2: 3179808, 3: 1024434, 4: 447984, 5: 94794, 8: 35136, 48: 186782}, 437248,
                  'kappa 0.001730418899766071163804 - ... (PR body)'),
    '3682 T360': ({1: 7980120, 2: 2345760, 3: 755730, 4: 330480, 5: 69930, 8: 25920, 48: 137790}, 322560,
                  'kappa 0.001730418899766071163804 - ... (PR body)'),
    '3815 T360': ({1: 7989840, 2: 2349000, 3: 749088, 4: 330480, 5: 69120, 8: 25920, 48: 137952}, 322560,
                  'kappa 0.001730578437266246818740 - ... (PR body; profile = #3945 old-PROFILE.json)'),
    '3945 T360': ({1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 48: 137952}, 322560,
                  'root [0.001730663236534549309146, ...147] (RESULT.json)'),
}
for name, (H, W, stated) in WORDS.items():
    if H is None:
        print(f'{name:10s} profile not available from the PR body; stated {stated}'); continue
    mass = sum(r * n for r, n in H.items())
    b = price_b(H, 73, W)
    b2 = root(H, {73: W})
    print(f'{name:10s} W={W} calls/children={sum(H.values())} mass={mass} deficit={73*W-mass} ({(73*W-mass)/(W//896)} per replica)')
    print(f'           hub-engine b (1e-24 grid, certified) = {dec(b, 25)}')
    print(f'           mixed root (Decimal)                 = {b2:.22e}')
    print(f'           stated: {stated}')
