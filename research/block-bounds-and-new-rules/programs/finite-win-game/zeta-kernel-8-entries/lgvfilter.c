/* lgvfilter.c -- exact minor test (Lindstrom-Gessel-Viennot) for call structures.
 * A word x_i <- x_i + c_t x_j (t = 1..G) is a weighted acyclic network (layers = times; stay edges weight 1,
 * jump edges j->i weight c_t).  By LGV (valid for any acyclic network), every minor Pi[R,C] equals
 *     sum over VERTEX-DISJOINT path systems from inputs C to final registers R of  sign(sigma) * weight.
 * Hence: if exactly ONE disjoint system exists, the minor is +- a single monomial, NONZERO for all c != 0.
 * A correct word has Pi = P D Z, so Pi[R,C] = 0  iff  Z[pi(R),C] is singular.  Necessary conditions:
 *     Z[pi(R),C] singular     ->  #systems(C,R) != 1
 *     Z[pi(R),C] nonsingular  ->  #systems(C,R) >= 1
 * (the 1x1 case is the 2-path rule for entries).  Reads leaves "pi ... | seq i<j ..." on stdin and prints the
 * survivors; counts on stderr.  Counting by DP over register subsets (counts capped at 2).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int N = 8, FULL = 255, ZR[8], RK[256][256];

static int rank_rows(int rowmask, int colmask) {
    long long M[8][8]; int m = 0, n = 0, cols[8];
    for (int u = 0; u < N; u++) if (colmask >> u & 1) cols[n++] = u;
    for (int s = 0; s < N; s++) if (rowmask >> s & 1) { for (int k = 0; k < n; k++) M[m][k] = (ZR[s] >> cols[k]) & 1; m++; }
    if (!m || !n) return 0;
    int r = 0; long long prev = 1;
    for (int c = 0; c < n && r < m; c++) {
        int p = -1; for (int i = r; i < m; i++) if (M[i][c]) { p = i; break; }
        if (p < 0) continue;
        if (p != r) for (int k = 0; k < n; k++) { long long t = M[p][k]; M[p][k] = M[r][k]; M[r][k] = t; }
        for (int i = r + 1; i < m; i++) {
            for (int k = c + 1; k < n; k++) M[i][k] = (M[r][c] * M[i][k] - M[i][c] * M[r][k]) / prev;
            M[i][c] = 0;
        }
        prev = M[r][c]; r++;
    }
    return r;
}

int main(int argc, char **argv) {
    if (argc > 1) { int b = atoi(argv[1]); N = 1 << b; FULL = (1 << N) - 1; }
    for (int s = 0; s < N; s++) { ZR[s] = 0; for (int u = 0; u < N; u++) if ((s & ~u) == 0) ZR[s] |= 1 << u; }
    for (int R = 0; R <= FULL; R++) for (int C = 0; C <= FULL; C++) RK[R][C] = rank_rows(R, C);
    char line[2048]; long long nin = 0, nout = 0, fail_sing = 0, fail_nons = 0;
    static unsigned char cnt[256];
    while (fgets(line, sizeof line, stdin)) {
        if (strncmp(line, "pi", 2)) continue;
        nin++;
        int pi[8], k = 0; char *p = line + 2;
        while (k < N) { char *e; long x = strtol(p, &e, 10); if (e == p) break; pi[k++] = (int)x; p = e; }
        int ci[64], cj[64], g = 0;
        char *q = strchr(line, '|');
        while (q && *q) { if (*q >= '0' && *q <= '9' && q[1] == '<') { ci[g] = *q - '0'; cj[g] = q[2] - '0'; g++; q += 3; } else q++; }
        int RM[256];
        for (int X = 0; X <= FULL; X++) { RM[X] = 0; for (int x = 0; x < N; x++) if (X >> x & 1) RM[X] |= 1 << pi[x]; }
        int ok = 1;
        for (int C = 1; C <= FULL && ok; C++) {
            memset(cnt, 0, sizeof cnt);
            cnt[C] = 1;
            for (int t = 0; t < g; t++) {
                int i = ci[t], j = cj[t];
                for (int S = 0; S <= FULL; S++) {
                    if ((S >> i & 1) && !(S >> j & 1)) {
                        int v = cnt[S] + cnt[S ^ (1 << i) ^ (1 << j)];
                        cnt[S] = v > 2 ? 2 : v;
                    }
                }
            }
            int sz = __builtin_popcount(C);
            for (int R = 1; R <= FULL; R++) {
                if (__builtin_popcount(R) != sz) continue;
                int sing = RK[RM[R]][C] < sz;
                if (sing && cnt[R] == 1) { ok = 0; fail_sing++; break; }
                if (!sing && cnt[R] == 0) { ok = 0; fail_nons++; break; }
            }
        }
        if (ok) { nout++; fputs(line, stdout); }
    }
    fprintf(stderr, "lgvfilter: in %lld, survivors %lld (killed: one-system singular minor %lld, no-system nonsingular minor %lld)\n",
            nin, nout, fail_sing, fail_nons);
    return 0;
}
