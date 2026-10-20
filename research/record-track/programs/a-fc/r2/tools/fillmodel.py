# fc-reopt round 2: the literal model of a #3945-line fill, read exactly as #3945's owner/verify.py reads it
# (windows, families, widths old/new, bins). Our code; the window/family/width logic is copied from verify.py.
import json, gzip
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np


def read(path):
    path = Path(path); raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)


def dimension(rows):
    piv = {}
    for v in rows:
        while v:
            b = v.bit_length()
            if b not in piv: piv[b] = v; break
            v ^= piv[b]
    return len(piv)


def canonical_widths(word, origin):
    full = [1 << i for i in range(9)]
    return [dimension(word['helper_final'].get(str(r), full)) - dimension(word['start_extra'].get(str(r), [])) for r in origin]


def load(owner_dir):
    """owner_dir = #3945's research/staged544-terminal-gap-fill/owner (inputs/ with old and new words)."""
    R = Path(owner_dir); OLD = R / 'inputs'; T = 360; N = 9 * T
    order = read(OLD / 'ORDER.json'); dispatch = read(OLD / 'DISPATCH.json')
    clock = {tuple(call): t for t, call in enumerate(order)}
    provider = {(s, r): tag for s, r, tag in dispatch}
    template = read(OLD / 'TEMPLATE.json.gz'); plan = read(OLD / 'PLAN.json.gz')
    children = {row['child']: row['rank'] for row in template['idle_children']}; K = max(children) + 1
    pair = {row['stage']: row['pair'] for row in template['recipe'] if row['kind'] == 'local'}
    windows = {}
    for owner in range(T):
        touched = {}; tick = None
        for event in template['execution_schedule']:
            kind = event['kind']
            if kind in ('bank_padding_cover_sweep', 'translate_every_live_family'): continue
            if kind == 'local_cover_sweep': tick = clock[event['stage'], owner]; touched[pair[event['stage']]] = tick
            elif kind == 'bridge': touched[event['src']] = tick; touched[event['dst']] = tick
            elif kind == 'idle_cover_sweep':
                for child in event['child_indices']:
                    if child in children:
                        assert (owner, child) not in windows
                        windows[owner, child] = (children[child], touched.get(event['pair'], -1), tick)
    family = {c: 2 * plan['pairs'][str(c)] + plan['sides'][str(c)] for c in children}
    widths = {}; oldwidths = {}; local = {}
    for tag, pid in [('base', 'gauge44'), ('extra', 'extra')]:
        word = read(R / f'inputs/WORD-{tag}.json.gz'); oldword = read(OLD / f'old-WORD-{tag}.json.gz')
        origin = read(OLD / f'old-NATIVE-{tag}.json.gz')['helper_origin']
        widths[pid] = canonical_widths(word, origin); oldwidths[pid] = canonical_widths(oldword, origin)
        local[pid] = {int(r): n for r, n in word['report']['local_profile'].items()}
    # per call: widths of its 408 roles (old, new)
    W_old = np.array([oldwidths[provider[tuple(order[t])]] for t in range(N)], dtype=np.int16)
    W_new = np.array([widths[provider[tuple(order[t])]] for t in range(N)], dtype=np.int16)
    return dict(T=T, N=N, order=order, provider=provider, children=children, K=K, windows=windows, family=family,
                W_old=W_old, W_new=W_new, local=local)


def bin_index(M, o, c, p):
    return (o * M['K'] + c) * 112 + p
