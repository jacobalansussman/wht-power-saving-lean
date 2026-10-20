/* box_enum2.c (ex2-families): ALL move multisets for a commutative-slice target (m, W, M) that pass the deletion test,
   by depth-first search over the direction counts c_1..c_n (n = 2^m - 1) inside the LP box [lo_z, hi_z] (read from
   stdin), with pruning: for every odd form Q, (weight so far) + (largest possible remaining weight) >= W.
   Prints each deletion-passing count vector (one line: c_1 .. c_n).  Input: m W M, then n pairs lo hi.
   Odd forms: binary quadratic forms Q(z) = sum_{j<=k} q_jk z_j z_k with sum_j q_jj = 1 (mod 2). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int m, W, M, n, P, nodd;
static int lo[64], hi[64], c[64];
static unsigned char *val;   /* nodd x n */
static int *w;               /* nodd */
static int *remmax;          /* (n+1) x nodd: max remaining weight from z onward given hi (ignores sum) */
static long long npass = 0, nodes = 0;
static void rec(int z, int rem) {
  nodes++;
  if (z == n) {
    if (rem) return;
    for (int q = 0; q < nodd; q++) if (w[q] < W) return;
    npass++;
    for (int i = 0; i < n; i++) printf("%d ", c[i]);
    printf("\n"); fflush(stdout);
    return;
  }
  /* sum feasibility */
  int slo = 0, shi = 0; for (int i = z; i < n; i++) { slo += lo[i]; shi += hi[i]; }
  if (rem < slo || rem > shi) return;
  for (int x = lo[z]; x <= hi[z]; x++) {
    if (x > rem) break;
    c[z] = x;
    int ok = 1;
    for (int q = 0; q < nodd; q++) {
      w[q] += x * val[q * n + z];
      int r = remmax[(z + 1) * nodd + q]; if (r > rem - x) r = rem - x;
      if (w[q] + r < W) ok = 0;
    }
    if (ok) rec(z + 1, rem - x);
    for (int q = 0; q < nodd; q++) w[q] -= x * val[q * n + z];
  }
  c[z] = 0;
}
int main(void) {
  if (scanf("%d %d %d", &m, &W, &M) != 3) return 1;
  n = (1 << m) - 1;
  for (int i = 0; i < n; i++) if (scanf("%d %d", &lo[i], &hi[i]) != 2) return 1;
  int pj[64], pk[64]; P = 0;
  for (int j = 0; j < m; j++) for (int k = j; k < m; k++) { pj[P] = j; pk[P] = k; P++; }
  long nf = 1L << P;
  val = malloc((size_t)nf * n); nodd = 0;
  for (long f = 0; f < nf; f++) {
    int tr = 0; for (int b = 0; b < P; b++) if ((f >> b & 1) && pj[b] == pk[b]) tr ^= 1;
    if (!tr) continue;
    for (int z = 1; z <= n; z++) {
      int v = 0; for (int b = 0; b < P; b++) if ((f >> b & 1) && (z >> pj[b] & 1) && (z >> pk[b] & 1)) v ^= 1;
      val[(size_t)nodd * n + (z - 1)] = v;
    }
    nodd++;
  }
  w = calloc(nodd, sizeof(int));
  remmax = calloc((size_t)(n + 1) * nodd, sizeof(int));
  for (int z = n - 1; z >= 0; z--) for (int q = 0; q < nodd; q++)
    remmax[z * nodd + q] = remmax[(z + 1) * nodd + q] + hi[z] * val[q * n + z];
  rec(0, M);
  fprintf(stderr, "m=%d W=%d M=%d: odd forms %d, nodes %lld, deletion-passing multisets %lld\n", m, W, M, nodd, nodes, npass);
  return 0;
}
