# The stop-band table of the report: every coefficient set on the re-centred prototype Bs, 2^21-point response.
import os, sys, json, numpy as np
os.environ['PROTO'] = 'Bs'; sys.argv = sys.argv[:1]; import gq
rows = [  # label, scheme or None (plain), s, file or ('round', w), DSP per signal-lane, extractions per pair-lane
 ('baseline: 1 product/DSP, 18-bit (as built)', None, None, ('round', 18), 14, 0),
 ('baseline: 1 product/DSP, 12-bit', None, None, ('round', 12), 14, 0),
 ('WP486-style per-tap packing, s=19 (12-bit)', 'pt', 19, 'c3_Bs_pt_s19.npy', 7, 14),
 ('ap-4 whole branch, s=19', 'whole', 19, 'c3_Bs_whole_s19.npy', 7, 1),
 ('ap-4 groups of 7, s=19', 'g7', 19, 'c3_Bs_g7_s19.npy', 7, 2),
 ('as-4 graded c2s, s=19 (DSP48E2)', 'c2s', 19, 'c3_Bs_c2s_s19.npy', 7, 2),
 ('as-4 graded g3, s=19 (DSP48E2)', 'g3', 19, 'c3_Bs_g3_s19.npy', 7, 3),
 ('WP486-style per-tap packing, s=17 (DSP48E1)', 'pt', 17, 'c3_Bs_pt_s17.npy', 7, 14),
 ('ap-4 whole branch, s=17 (DSP48E1)', 'whole', 17, 'c3_Bs_whole_s17.npy', 7, 1),
 ('ap-4 groups of 7, s=17 (DSP48E1)', 'g7', 17, 'c3_Bs_g7_s17.npy', 7, 2),
 ('as-4 graded g4, s=17 (DSP48E1)', 'g4', 17, os.environ.get('G4', 'c3_Bs_g4_s17.npy'), 7, 4),
 ('as-4 graded g5, s=17 (DSP48E1)', 'g5', 17, 'final_Bs_g5_s17.npy', 7, 5),
 ('as-4 u4: centre pair unpacked + graded, s=17 (DSP48E1)', 'u4', 17, 'c3_Bs_u4_s17.npy', 8, 3),
]
print("| design (prototype Bs, linear phase) | DSP per signal-lane | extractions per pair-lane | ripple ±dB (need <= 0.2) | alias edge [0.685, 1.5) ch dB (need <= -60) | worst beyond 1.5 ch dB (target <= -80) | median beyond 1.5 ch | share of stop-band <= -80 dB | legal | file |")
print("|---|---|---|---|---|---|---|---|---|---|")
for lab, sc, s, f, dsp, ex in rows:
    if isinstance(f, tuple):
        c = np.round(gq.h / np.abs(gq.h).max() * ((1 << (f[1] - 1)) - 1)).astype(np.int64); ok = '-'; fn = f'rounded {f[1]}-bit'
    else:
        if os.environ.get('C2') and os.path.exists(f.replace('.npy', '_c2.npy')): f = f.replace('.npy', '_c2.npy')
        if not os.path.exists(f): print(f"| {lab} | missing {f} |"); continue
        c = np.load(f); sh, gid = gq.setup(sc); ok = gq.legal(gq.raw_of(c, sh), gid, s); fn = f
    m = gq.metrics(c)
    print(f"| {lab} | {dsp} | {ex} | {m['ripple']:.3f} | {m['stop1']:.2f} | {m['stop2']:.2f} | {m['med2']:.1f} | {100*m['frac80']:.1f} % | {ok} | `{fn}` |")
