/* tdfs.c -- exhaustive structural search for in-place words computing Z = T^(x)b, T=[[1,1],[0,1]].
 *
 * Model (OpenAI finite-win game for A = T): N = 2^b registers, register c starts with x_c.  A call (i <- j)
 * replaces x_i by a*x_i + b*x_j, a,b != 0; monomials free.  WLOG (push monomials to the ends, conjugate) every
 * call is x_i <- x_i + c x_j (c != 0), the right monomial is the identity (register c starts with input c), and the
 * end state is P*D*Z: register i ends with d_i * row pi(i) of Z, row s = sum_{U >= s} x_U.
 *
 * For a fixed pi (read from stdin, one per line), enumerate call sequences of length G in "adjacent normal form"
 * (adjacent commuting calls in increasing index; commuting = neither writes what the other reads) and prune with
 * NECESSARY conditions only:
 *   (cut)  for every register set X and time t: future calls i<-j with i in X, j notin X  >=
 *          rank Z[pi(X), complement(K_X(t))], K_X(t) = union of structural knowledge sets of X.
 *          [Grassmann distance: only such calls move span(rows of X), each by one step]  Checked through the best set-partition (each call enters at most
 *          one block of a partition), computed by DP over subsets.
 *   (zero) an entry (i,u) with u in K_i but u notin supp(row pi(i)) needs >= 2 paths at the end (else it is a
 *          single monomial, nonzero); while it has 1 path, register i needs another call: counted in the singleton.
 * Leaf (t = G): K_i >= supp row pi(i) for all i; every extra entry has >= 2 paths; and at least one call had
 *          overlapping knowledge (K_i & K_j != 0) -- required only with flag NEEDCANCEL (a word with no
 *          overlapping call has single paths everywhere, hence is cancellation-free).
 * Output: every surviving leaf (pi + sequence) to stdout (capped), counts to stderr.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

static int B, N, G, FULL, NEEDCANCEL, LGV;
static long long MAXPRINT;
static int pi_[8], Tm[8];
static int RK[256][256];
static int seqi[32], seqj[32];
static long long nodes, leaves, printed, MAXNODES; static int aborted;
static int ZR[8];   /* row s of Z as a column mask */

static int rank_rows(int rowmask, int colmask) {
    /* exact rank over Q of the 0/1 submatrix Z[rows][cols] via fraction-free (Bareiss) elimination */
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

static int needX[256], best[256], RMX[256]; static double YD[256], YBOUND; static int HAVEY;

static int partition_bound(const int *K, const int *bad) {
    static int km[256];
    km[0] = 0;
    for (int X = 1; X <= FULL; X++) {
        int low = __builtin_ctz(X);
        km[X] = km[X & (X - 1)] | K[low];
        needX[X] = RK[RMX[X]][FULL & ~km[X]];
    }
    for (int i = 0; i < N; i++) if (bad[i] && needX[1 << i] < 1) needX[1 << i] = 1;
    if (HAVEY) { double sy = 0; for (int X = 1; X < FULL; X++) sy += YD[X] * needX[X]; YBOUND = sy; }
    best[0] = 0;
    for (int S = 1; S <= FULL; S++) {
        int low = S & -S, rest = S ^ low, bb = 0;
        for (int sub = rest;; sub = (sub - 1) & rest) {
            int X = sub | low, v = needX[X] + best[S ^ X];
            if (v > bb) bb = v;
            if (!sub) break;
        }
        best[S] = bb;
    }
    return best[FULL];
}


/* DUAL (inverse-transpose) bound: Pi^{-T} evolves by the same calls with arcs reversed (call i<-j makes
 * K*_j |= K*_i), from I to P D^{-1} Z^{-T}; register r ends with column pi(r) of Z^{-1} (support {s <= pi(r)}).
 * Grassmann lemma for that process: future calls i<-j with j in X, i notin X  >=  rank Z[ complement(K*_X), pi(X) ]
 * (sign scalings do not change ranks).  Each call leaves exactly one block of a partition. */
static int needD[256], bestD[256];
static int dual_partition_bound(const int *KD, const int *badD) {
    static int km[256];
    km[0] = 0;
    for (int X = 1; X <= FULL; X++) {
        int low = __builtin_ctz(X);
        km[X] = km[X & (X - 1)] | KD[low];
        needD[X] = RK[FULL & ~km[X]][RMX[X]];
    }
    for (int i = 0; i < N; i++) if (badD[i] && needD[1 << i] < 1) needD[1 << i] = 1;
    bestD[0] = 0;
    for (int S = 1; S <= FULL; S++) {
        int low = S & -S, rest = S ^ low, bb = 0;
        for (int sub = rest;; sub = (sub - 1) & rest) {
            int X = sub | low, v = needD[X] + bestD[S ^ X];
            if (v > bb) bb = v;
            if (!sub) break;
        }
        bestD[S] = bb;
    }
    return bestD[FULL];
}
static int DUAL, Dm[8];

static int CORD[256], NCORD;
static long long lgv_killed;
static int lgv_ok(void) {
    static unsigned char cnt[256];
    for (int a = 0; a < NCORD; a++) {
        int C = CORD[a], sz = __builtin_popcount(C);
        memset(cnt, 0, sizeof cnt);
        cnt[C] = 1;
        for (int t = 0; t < G; t++) {
            int i = seqi[t], j = seqj[t];
            for (int S = 0; S <= FULL; S++)
                if ((S >> i & 1) && !(S >> j & 1)) { int v = cnt[S] + cnt[S ^ (1 << i) ^ (1 << j)]; cnt[S] = v > 2 ? 2 : v; }
        }
        for (int R = 1; R <= FULL; R++) {
            if (__builtin_popcount(R) != sz) continue;
            int sing = RK[RMX[R]][C] < sz;
            if (sing ? cnt[R] == 1 : cnt[R] == 0) return 0;
        }
    }
    return 1;
}

static int KDs[8], M2Ds[8];
static void dfs(int t, int *K, int *M2, int overlap, int pi_prev, int pj_prev) {
    int *KD = KDs, *M2D = M2Ds;
    nodes++;
    if (MAXNODES && nodes > MAXNODES) { aborted = 1; return; }
    int bad[8];
    for (int i = 0; i < N; i++) bad[i] = ((K[i] & ~Tm[i] & ~M2[i]) != 0);
    int r = G - t;
    if (t == G) {
        for (int i = 0; i < N; i++) { if ((K[i] & Tm[i]) != Tm[i]) return; if (bad[i]) return; }
        if (NEEDCANCEL && !overlap) return;
        leaves++;
        if (LGV && !lgv_ok()) { lgv_killed++; return; }
        if (printed < MAXPRINT) {
            printed++;
            printf("pi");
            for (int i = 0; i < N; i++) printf(" %d", pi_[i]);
            printf(" | seq");
            for (int s = 0; s < G; s++) printf(" %d<%d", seqi[s], seqj[s]);
            printf("\n");
        }
        return;
    }
    int badD[8];
    if (DUAL) for (int i = 0; i < N; i++) badD[i] = ((KD[i] & ~Dm[i] & ~M2D[i]) != 0);
    if (t == G && DUAL) { for (int i = 0; i < N; i++) if ((KD[i] & Dm[i]) != Dm[i] || badD[i]) return; }
    { int s1 = 0; for (int i = 0; i < N; i++) s1 += (bad[i] || (K[i] & Tm[i]) != Tm[i]); if (s1 > r) return; }
    if (DUAL) { int s2 = 0; for (int i = 0; i < N; i++) s2 += (badD[i] || (KD[i] & Dm[i]) != Dm[i]); if (s2 > r) return; }
    if (partition_bound(K, bad) > r) return;
    if (HAVEY && YBOUND > r + 1e-6) return;
    if (DUAL && dual_partition_bound(KD, badD) > r) return;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        if (i == j) continue;
        if (t > 0) {
            int commute = (i != pj_prev) && (pi_prev != j);
            if (commute && (i * N + j) <= (pi_prev * N + pj_prev)) continue;
        }
        int oK = K[i], oM = M2[i], oKD = KD[j], oMD = M2D[j];
        int nk = K[i] | K[j];
        int nm = M2[i] | M2[j] | (K[i] & K[j]);
        K[i] = nk; M2[i] = nm;
        M2D[j] = M2D[j] | M2D[i] | (KD[j] & KD[i]); KD[j] = KD[j] | KD[i];
        seqi[t] = i; seqj[t] = j;
        dfs(t + 1, K, M2, overlap | ((oK & K[j]) != 0), i, j);
        K[i] = oK; M2[i] = oM; KD[j] = oKD; M2D[j] = oMD;
    }
}

/* walk mode: check that a GIVEN word (pi | seq) is never pruned; prints the first failing step, if any */
static int walk(const char *line) {
    int v[8], k = 0; const char *p = line;
    while (*p && *p != '|') { if (*p >= '0' && *p <= '9') v[k++] = *p - '0'; p++; }
    for (int i = 0; i < N; i++) { pi_[i] = v[i]; Tm[i] = ZR[v[i]]; }
    for (int X = 0; X <= FULL; X++) { RMX[X] = 0; for (int x = 0; x < N; x++) if (X >> x & 1) RMX[X] |= 1 << pi_[x]; }
    int ci[64], cj[64], g = 0;
    while (*p) { if (*p >= '0' && *p <= '9' && p[1] == '<') { ci[g] = *p - '0'; cj[g] = p[2] - '0'; g++; p += 3; } else p++; }
    int K[8], M2[8], overlap = 0;
    for (int i = 0; i < N; i++) { K[i] = 1 << i; M2[i] = 0; }
    for (int t = 0; t <= g; t++) {
        int bad[8], r = g - t;
        for (int i = 0; i < N; i++) bad[i] = ((K[i] & ~Tm[i] & ~M2[i]) != 0);
        if (t < g) {
            int pb = partition_bound(K, bad);
            if (pb > r) { printf("walk: PRUNED at t=%d (partition bound %d > remaining %d)\n", t, pb, r); return 1; }
            int i = ci[t], j = cj[t];
            overlap |= (K[i] & K[j]) != 0;
            M2[i] = M2[i] | M2[j] | (K[i] & K[j]); K[i] |= K[j];
        } else {
            for (int i = 0; i < N; i++) if ((K[i] & Tm[i]) != Tm[i] || bad[i]) { printf("walk: leaf condition FAILS at register %d\n", i); return 1; }
            printf("walk: word of %d calls passes every pruning test (overlap=%d)\n", g, overlap);
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: tdfs b G [needcancel=1] [maxprint=100000] < pis\n"); return 1; }
    B = atoi(argv[1]); G = atoi(argv[2]);
    NEEDCANCEL = argc > 3 ? atoi(argv[3]) : 1;
    MAXPRINT = argc > 4 ? atoll(argv[4]) : 100000;
    MAXNODES = argc > 5 ? atoll(argv[5]) : 0;
    LGV = argc > 6 ? atoi(argv[6]) : 0;
    DUAL = argc > 7 ? atoi(argv[7]) : 0;
    N = 1 << B; FULL = (1 << N) - 1;
    for (int s = 0; s < N; s++) { ZR[s] = 0; for (int u = 0; u < N; u++) if ((s & ~u) == 0) ZR[s] |= 1 << u; }
    for (int R = 0; R <= FULL; R++) for (int C = 0; C <= FULL; C++) RK[R][C] = rank_rows(R, C);
    NCORD = 0; for (int pc = 1; pc <= N; pc++) for (int C = 1; C <= FULL; C++) if (__builtin_popcount(C) == pc) CORD[NCORD++] = C;
    char line[16384]; long long totleaves = 0; int npi = 0;
    if (G == 0) { while (fgets(line, sizeof line, stdin)) walk(line); return 0; }
    while (fgets(line, sizeof line, stdin)) {
        int v[8], k = 0; char *p = line;
        while (k < N) { char *e; long x = strtol(p, &e, 10); if (e == p) break; v[k++] = (int)x; p = e; }
        if (k < N) continue;
        HAVEY = 0;
        { char *yy = strstr(line, "| y");
          if (yy) { yy += 3; HAVEY = 1; for (int X = 1; X < FULL; X++) { char *e; YD[X] = strtod(yy, &e); yy = e; } } }
        for (int i = 0; i < N; i++) { pi_[i] = v[i]; Tm[i] = ZR[v[i]]; }
        for (int X = 0; X <= FULL; X++) { RMX[X] = 0; for (int x = 0; x < N; x++) if (X >> x & 1) RMX[X] |= 1 << pi_[x]; }
        int K[8], M2[8];
        for (int i = 0; i < N; i++) { K[i] = 1 << i; M2[i] = 0; KDs[i] = 1 << i; M2Ds[i] = 0; }
        for (int i = 0; i < N; i++) { Dm[i] = 0; for (int s2 = 0; s2 < N; s2++) if ((s2 & ~pi_[i]) == 0) Dm[i] |= 1 << s2; }
        nodes = 0; leaves = 0; aborted = 0; printed = 0; lgv_killed = 0;
        clock_t c0 = clock();
        dfs(0, K, M2, 0, -1, -1);
        npi++; totleaves += leaves - lgv_killed;
        fprintf(stderr, "pi");
        for (int i = 0; i < N; i++) fprintf(stderr, " %d", pi_[i]);
        fprintf(stderr, " : nodes %lld leaves %lld lgv-killed %lld survivors %lld %s %.1fs\n", nodes, leaves, lgv_killed, leaves - lgv_killed, aborted ? "ABORTED(node cap)" : "complete", (double)(clock() - c0) / CLOCKS_PER_SEC);
        fflush(stdout); fflush(stderr);
    }
    fprintf(stderr, "TOTAL pis %d leaves %lld\n", npi, totleaves);
    return 0;
}
