"""Z7 per-character screens (t2-block3's own: hook 80 m_k >= 224 a + 32 b; Theorem D trivial m_0 >= a+b+2;
non-trivial m_k >= a+2 (a>0); depth m_k >= 3 when any piece is present) in two frameworks:
 EXACT: every step mapped exactly by tau (t2-block3): all m_k = j, M = 7j.
 GAUGE: register values mapped up to free letters; a move on a register whose content is one character's
        (right-eigen under c', trivial left gauge) is ONE move acting in that character only, so the m_k are independent, M = sum m_k."""
import math
def lb(a, b, k):
    h = math.ceil((224 * a + 32 * b) / 80)
    if k == 0: return max(h, a + b + 2)
    return max(h, (a + 2) if a else 0, 3 if (a + b) else 0)
def best(W, mode):
    res = None
    for a in range(W // 7 + 1):
        b = W - 7 * a
        m0, m1 = lb(a, b, 0), lb(a, b, 1)
        M = 7 * max(m0, m1) if mode == "EXACT" else m0 + 6 * m1
        if res is None or M < res[0]: res = (M, a, b, m0, m1)
    return res
print("W  plain  EXACT(minM,a,b)  GAUGE(minM,a,b,m0,m_k)   targets W in 3W-1 / 3W-2 pass?")
for W in list(range(2, 41)):
    e, g = best(W, "EXACT"), best(W, "GAUGE")
    tag = []
    for M in (3 * W - 1, 3 * W - 2):
        pe = any(7 * j == M2 for j in range(1, 200) for M2 in [M] if 7 * j >= e[0] and 7 * j <= M)
        pe = (e[0] <= M)   # EXACT survivors are M' = 7j with e[0] <= M' <= M
        pg = g[0] <= M
        tag.append("%d in %d: exact %s gauge %s" % (W, M, "pass" if pe else "-", "pass" if pg else "-"))
    if g[0] <= 3 * W - 1 or e[0] <= 3 * W - 1:
        print(W, 3 * W, e[:3], g, " | ".join(tag))
for (W, M) in [(5, 14), (7, 20), (10, 28), (10, 29), (12, 35), (15, 44), (21, 62), (20, 59), (20, 58), (15, 42), (6, 17)]:
    e, g = best(W, "EXACT"), best(W, "GAUGE")
    ex = any(7 * j <= M and 7 * j >= e[0] for j in range(200))
    print("target %d in %d: EXACT %s   GAUGE %s (gauge min M %d at a=%d b=%d m0=%d mk=%d)" % (W, M, "survives" if ex else "excluded", "survives" if g[0] <= M else "excluded", *g))
