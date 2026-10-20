"""xw3-extract: the extracted collapse vectors (fu_<tag>.npz, from fu.py on the Colabs), read back.
For every flag: the f=3,4 collapse space X (columns = vectors over the kept top words); its f=5 SURVIVOR part Xs
(when any); and for the space orthogonal to the survivors (the part that fails at f=5) the top words that carry
weight, with the letters they use.  Label of a top word: (mix before move 1 = '1', mix before move 2, ..., readout);
the word is  readout . M_n . mix . ... . M_2 . mix . M_1 applied to the input v.
"""
import sys, os, collections
import numpy as np
here = os.environ.get("FU_DIR", os.path.dirname(os.path.abspath(__file__)))  # folder with fu_<tag>.npz
SPECIAL = ("lX", "lS", "fS", "fX", "B", "Bi", "A", "Ai", "a0", "a0i", "E01", "Pr", "T2", "T0", "CS", "A.T2.Ai", "Ai.CS01.A")


def letters(lab):
    out = set()
    for w in lab:
        for x in w.split("."):
            for s in SPECIAL:
                if x == s or (s in ("lX", "lS", "fS", "fX", "Pr", "CS") and x.startswith(s)): out.add(x)
    return out


for tag in sys.argv[1:] or ["B1", "B2", "B3", "C1", "C2", "C3", "C4s601"]:
    d = np.load(os.path.join(here, f"fu_{tag}.npz"))
    keep, X, Xs = d["keep"], d["X"], d["Xs"]
    lab = [tuple(s.split("|")) for s in d["lab"]]
    nt = X.shape[1]; ns = Xs.shape[1] if Xs.ndim == 2 else 0
    print(f"=== {tag}: {nt} collapse vectors at f=3,4 over {keep.size} kept top words; {ns} survive at f=5")
    # orthonormal basis of span X; split off the survivor span
    Qx, _ = np.linalg.qr(X)
    if ns:
        Qs, _ = np.linalg.qr(Xs)
        F = Qx - Qs @ (Qs.conj().T @ Qx)
        U, s, _ = np.linalg.svd(F, full_matrices=False); Qf = U[:, s > 1e-8]
    else:
        Qf = Qx
    weight = (np.abs(Qf) ** 2).sum(axis=1)          # weight of each top word in the failing part
    order = np.argsort(-weight)
    tot = weight.sum()
    cum = np.cumsum(weight[order]) / tot
    n90 = int(np.searchsorted(cum, 0.9) + 1)
    print(f"  failing part: dim {Qf.shape[1]}; {n90} top words carry 90% of its weight; heaviest:")
    for i in order[:12]:
        print(f"    {lab[keep[i]]}  weight {weight[i] / tot:.3f}")
    cnt = collections.Counter()
    for i in order[:n90]:
        for l in letters(lab[keep[i]]): cnt[l] += weight[i] / tot
    print(f"  letters in the words carrying 90% (share of weight): "
          f"{', '.join(f'{k} {v:.2f}' for k, v in cnt.most_common(14))}")
    if ns:
        ws = (np.abs(np.linalg.qr(Xs)[0]) ** 2).sum(axis=1); o = np.argsort(-ws)
        print(f"  survivor part: dim {ns}; heaviest words: {[lab[keep[i]] for i in o[:10]]}")
