/* 23 Oct 2026.  GENERATOR (untrusted) of nonsingular-minor certificates for the hook demand of a PARTIAL product of
 * letters, g_L = prod over b in L of dir(e_b), L a set of address bits (the "lab" mask; L = all bits gives K and the
 * matrix of gen_m0.c).  Same construction as gen_m0.c with Q replaced by
 *   Q_L[u][v] = (-i)^{popcount(u xor v)} if u xor v lies inside L, else 0      (g_L = ((1+i)/2)^|L| Q_L).
 *   rows: k-subsets U of {1..n-1}; cols: k-subsets T of {0..n-1} with 0 in T;
 *   M_alpha[U][T] = sum_s (-1)^{alpha.s} Q_L[0][s] det[ Q_L[u_a][t_b xor s] ]_{a,b}.
 * The demand of g_L on the hook N, f(g_L) = rank (1 - P_N)(g_L (x) Lambda^k g_L) P_N, is the sum over the 2^m
 * characters alpha of rank M_alpha.  Computes M_alpha mod P, eliminates (pivot per column, first usable row), prints
 * the rank and writes the pivot rows/cols (bitmasks) to OUTFILE.
 * Build: cc -O2 -DP=2147483629ULL -o gen_mj gen_mj.c      Run: ./gen_mj m k lab alpha OUTFILE
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifndef P
#define P 2147483629ULL
#endif
typedef uint64_t u64;
static u64 pw(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return r; }
static int popc(unsigned x) { return __builtin_popcount(x); }
static u64 pmi[4]; static unsigned lab;
static u64 q(unsigned u, unsigned v) { unsigned d = u ^ v; return (d & ~lab) ? 0 : pmi[popc(d) & 3]; }
int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: gen_mj m k lab alpha OUTFILE\n"); return 2; }
    int m = atoi(argv[1]), k = atoi(argv[2]); lab = (unsigned)atoi(argv[3]); unsigned alpha = (unsigned)atoi(argv[4]);
    int n = 1 << m; lab &= (unsigned)(n - 1);
    u64 I = 0;
    for (u64 g = 2; g < 1000; g++) if (pw(g, (P - 1) / 2) == P - 1) { I = pw(g, (P - 1) / 4); break; }
    if (I * I % P != P - 1) { fprintf(stderr, "no sqrt(-1)\n"); return 1; }
    pmi[0] = 1; pmi[1] = (P - I) % P; pmi[2] = P - 1; pmi[3] = I;          /* (-i)^0..3 */
    int NR = 0, NC = 0; unsigned *rows = malloc(sizeof(unsigned) << n), *cols = malloc(sizeof(unsigned) << n);
    for (unsigned mask = 0; mask < (1u << n); mask++) if (popc(mask) == k) {
        if (!(mask & 1)) rows[NR++] = mask; else cols[NC++] = mask;
    }
    fprintf(stderr, "m=%d k=%d lab=%u alpha=%u: %d rows, %d cols, P=%llu, I=%llu\n", m, k, lab, alpha, NR, NC,
            (unsigned long long)P, (unsigned long long)I);
    unsigned *smask = malloc(sizeof(unsigned) * NC * n); signed char *ssign = malloc(NC * n);
    for (int c = 0; c < NC; c++) {
        int t[32], L = 0; for (int x = 0; x < n; x++) if (cols[c] >> x & 1) t[L++] = x;
        for (int s = 0; s < n; s++) {
            int seq[32], inv = 0; unsigned mk = 0;
            for (int b = 0; b < k; b++) { seq[b] = t[b] ^ s; mk |= 1u << seq[b]; }
            for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) if (seq[a] > seq[b]) inv++;
            smask[c * n + s] = mk; ssign[c * n + s] = (inv & 1) ? -1 : 1;
        }
    }
    int *cnt = calloc(n + 1, sizeof(int)); unsigned **byp = malloc(sizeof(unsigned *) * (n + 1));
    for (unsigned mask = 0; mask < (1u << n); mask++) cnt[popc(mask)]++;
    for (int j = 0; j <= n; j++) { byp[j] = malloc(sizeof(unsigned) * (cnt[j] + 1)); cnt[j] = 0; }
    for (unsigned mask = 0; mask < (1u << n); mask++) { int j = popc(mask); byp[j][cnt[j]++] = mask; }
    u64 *D = malloc(sizeof(u64) << n);
    u64 *M = malloc(sizeof(u64) * (size_t)NR * NC);
    if (!M || !D) { fprintf(stderr, "alloc\n"); return 1; }
    for (int r = 0; r < NR; r++) {
        int u[32], L = 0; for (int x = 0; x < n; x++) if (rows[r] >> x & 1) u[L++] = x;
        /* D[W] = det Q_L[u_0..u_{j-1}, W] for |W| = j, Laplace along the last row */
        D[0] = 1;
        for (int j = 1; j <= k; j++) {
            for (int qi = 0; qi < cnt[j]; qi++) {
                unsigned W = byp[j][qi]; u64 acc = 0; int pos = 0;
                for (int w = 0; w < n; w++) if (W >> w & 1) {
                    u64 term = q((unsigned)u[j - 1], (unsigned)w) * D[W & ~(1u << w)] % P;
                    if (((j - 1) + pos) & 1) acc += P - term; else acc += term;
                    pos++;
                }
                D[W] = acc % P;
            }
        }
        u64 *row = M + (size_t)r * NC;
        for (int c = 0; c < NC; c++) {
            u64 acc = 0;
            for (int s = 0; s < n; s++) {
                u64 q0 = q(0, (unsigned)s); if (!q0) continue;
                u64 term = q0 * D[smask[c * n + s]] % P;
                int sg = ssign[c * n + s] * ((popc(alpha & s) & 1) ? -1 : 1);
                acc += (sg > 0) ? term : (P - term) % P;
            }
            row[c] = acc % P;
        }
    }
    fprintf(stderr, "matrix built\n");
    char *used = calloc(NR, 1); int *pr = malloc(sizeof(int) * NC), *pc = malloc(sizeof(int) * NC); int rank = 0;
    for (int c = 0; c < NC; c++) {
        int p0 = -1;
        for (int r = 0; r < NR; r++) if (!used[r] && M[(size_t)r * NC + c]) { p0 = r; break; }
        if (p0 < 0) continue;
        used[p0] = 1; pr[rank] = p0; pc[rank] = c; rank++;
        u64 *rp = M + (size_t)p0 * NC; u64 inv = pw(rp[c], P - 2);
        for (int r = 0; r < NR; r++) {
            if (used[r]) continue;
            u64 *rr = M + (size_t)r * NC; if (!rr[c]) continue;
            u64 f = (P - rr[c] * inv % P) % P;
            for (int j = c; j < NC; j++) rr[j] = (rr[j] + f * rp[j]) % P;
        }
        if ((rank & 511) == 0) fprintf(stderr, "  col %d rank %d\n", c, rank);
    }
    printf("m=%d k=%d lab=%u alpha=%u: rank mod %llu of M_alpha (%d x %d) = %d\n", m, k, lab, alpha,
           (unsigned long long)P, NR, NC, rank);
    FILE *fo = fopen(argv[5], "w");
    fprintf(fo, "m %d k %d lab %u alpha %u p %llu I %llu size %d\n", m, k, lab, alpha, (unsigned long long)P,
            (unsigned long long)I, rank);
    for (int i = 0; i < rank; i++) fprintf(fo, "%u %u\n", rows[pr[i]], cols[pc[i]]);
    fclose(fo);
    return 0;
}
