/* rcdp.c (rm2-exact, 15 Oct 2026): guillotine rank-cover bound for the multi-target system, fixed move counts.

For a form F with support D_F (a set of directions; moves = directions with multiplicity counts[z]):
   rank E[D_F, D_F] = s_F  (exact; Theorem 10 with a target per array, rm-theory 8.5).
For R, C subsets of D_F:  rank E[R, C] <= min(|R|, |C|, s_G for every form G with R u C inside D_G), and
rank E[R0 u R1, C] <= rank E[R0, C] + rank E[R1, C] (same for columns).  The DP minimises this bound over guillotine
splits by linear functionals (so R, C range over D_F n affine flats).  If the minimum is < s_F: no program (exact).
Input  (binary, little endian int32): m nd nf ntest dmax ; counts[nd] ; mask[nf] ; s[nf] ; test[ntest]
Output (text): one line per tested form "F d s_F bound" and, for every violation, the certificate tree:
   "CERT F" then lines "P rowmask colmask bound how G"  (how: 0 = dimension, 1 = form G), terminated by "END".
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int m, nd, nf, ntest, dmax;
static int *counts, *smat, *tests; static uint32_t *mask;
static int d; static int dirs[32];
static int *best; static int *bestG;           /* over subsets of D_F (compressed) */
static int nsets; static uint32_t *sets; static int *wset; static int *split0, *split1;  /* split[i*nd + phi] */
static int *memo; static int *how;   /* memo[i*nsets+j] */
#define INF 1000000000
#define HMASK ((1<<22)-1)
static int *htab; static uint32_t *hkey;

static uint32_t compress(uint32_t dm) {   /* direction mask -> D_F bit mask */
  uint32_t r = 0; for (int k = 0; k < d; k++) if (dm >> (dirs[k]-1) & 1) r |= 1u << k; return r; }
static uint32_t expand(uint32_t cm) { uint32_t r = 0; for (int k = 0; k < d; k++) if (cm >> k & 1) r |= 1u << (dirs[k]-1); return r; }

static int find_or_add(uint32_t key) {
  uint32_t h = (key * 2654435761u) & HMASK;
  while (htab[h] >= 0) { if (hkey[h] == key) return htab[h]; h = (h + 1) & HMASK; }
  htab[h] = nsets; hkey[h] = key; sets[nsets] = key;
  int w = 0; for (int k = 0; k < d; k++) if (key >> k & 1) w += counts[dirs[k]-1];
  wset[nsets] = w; return nsets++;
}

static int cost(int i, int j) {
  if (i > j) { int t = i; i = j; j = t; }
  long idx = (long)i * nsets + j;
  if (memo[idx] >= 0) return memo[idx];
  int wi = wset[i], wj = wset[j];
  if (wi == 0 || wj == 0) { memo[idx] = 0; how[idx] = 0; return 0; }
  int c = wi < wj ? wi : wj; int hw = 0;
  uint32_t u = sets[i] | sets[j];
  if (best[u] < c) { c = best[u]; hw = 1; }
  if (c > 0) {
    for (int ph = 0; ph < nd; ph++) {
      int a = split0[i * nd + ph], b = split1[i * nd + ph];
      if (wset[a] && wset[b]) { int v = cost(a, j); if (v < c) { v += cost(b, j); if (v < c) { c = v; hw = 2 + 2 * ph; } } }
      a = split0[j * nd + ph]; b = split1[j * nd + ph];
      if (wset[a] && wset[b]) { int v = cost(i, a); if (v < c) { v += cost(i, b); if (v < c) { c = v; hw = 3 + 2 * ph; } } }
    }
  }
  memo[idx] = c; how[idx] = hw; return c;
}

static void cert(int ri, int rj) {  /* piece rows = set ri, columns = set rj */
  if (wset[ri] == 0 || wset[rj] == 0) return;
  int a = ri < rj ? ri : rj, b = ri < rj ? rj : ri; long idx = (long)a * nsets + b;
  int hw = how[idx], c = memo[idx];
  if (hw == 0) { printf("P %u %u %d 0 -1\n", expand(sets[ri]), expand(sets[rj]), c); return; }
  if (hw == 1) { printf("P %u %u %d 1 %d\n", expand(sets[ri]), expand(sets[rj]), c, bestG[sets[ri] | sets[rj]]); return; }
  int ph = (hw - 2) / 2; int s = ((hw - 2) % 2 == 0) ? a : b;
  int s0 = split0[s * nd + ph], s1 = split1[s * nd + ph];
  if (s == ri) { cert(s0, rj); cert(s1, rj); } else { cert(ri, s0); cert(ri, s1); }
}

int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "rb"); int hdr[5];
  if (!f || fread(hdr, 4, 5, f) != 5) { fprintf(stderr, "bad input\n"); return 1; }
  m = hdr[0]; nd = hdr[1]; nf = hdr[2]; ntest = hdr[3]; dmax = hdr[4];
  counts = malloc(4 * nd); mask = malloc(4 * nf); smat = malloc(4 * nf); tests = malloc(4 * (ntest + 1));
  if (fread(counts, 4, nd, f) != (size_t)nd || fread(mask, 4, nf, f) != (size_t)nf || fread(smat, 4, nf, f) != (size_t)nf
      || fread(tests, 4, ntest, f) != (size_t)ntest) { fprintf(stderr, "short input\n"); return 1; }
  fclose(f);
  best = malloc(sizeof(int) << dmax); bestG = malloc(sizeof(int) << dmax);
  int maxsets = 1 << 16;
  sets = malloc(4 * maxsets); wset = malloc(4 * maxsets);
  split0 = malloc(4L * maxsets * nd); split1 = malloc(4L * maxsets * nd);
  htab = malloc(4 * (HMASK + 1)); hkey = malloc(4 * (HMASK + 1));
  int nviol = 0;
  for (int t = 0; t < ntest; t++) {
    int F = tests[t]; d = 0;
    for (int z = 1; z <= nd; z++) if (mask[F] >> (z - 1) & 1) dirs[d++] = z;
    if (d > dmax || d == 0) { printf("F %d d %d skip\n", F, d); continue; }
    int sF = smat[F];
    uint32_t full = (d == 32) ? 0xffffffffu : ((1u << d) - 1);
    for (uint32_t u = 0; u <= full; u++) { best[u] = INF; bestG[u] = -1; if (u == full) break; }
    for (int G = 0; G < nf; G++) {
      if (G == F) continue;
      uint32_t u = compress(mask[G] & mask[F]);
      if (smat[G] < best[u]) { best[u] = smat[G]; bestG[u] = G; }
    }
    for (int b = 0; b < d; b++) for (uint32_t u = 0; u <= full; u++) {
      if (!(u >> b & 1)) { uint32_t v = u | (1u << b); if (best[v] < best[u]) { best[u] = best[v]; bestG[u] = bestG[v]; } }
      if (u == full) break; }
    /* the family of D_F n affine flats, closed under splitting by functionals */
    memset(htab, -1, 4 * (HMASK + 1)); nsets = 0;
    find_or_add(full); find_or_add(0);
    for (int i = 0; i < nsets; i++) {
      for (int ph = 1; ph <= nd; ph++) {
        uint32_t s0 = 0, s1 = 0;
        for (int k = 0; k < d; k++) if (sets[i] >> k & 1) { if (__builtin_popcount(ph & dirs[k]) & 1) s1 |= 1u << k; else s0 |= 1u << k; }
        int a = find_or_add(s0), b = find_or_add(s1);
        if (nsets >= maxsets) { fprintf(stderr, "too many sets\n"); return 2; }
        split0[i * nd + ph - 1] = a; split1[i * nd + ph - 1] = b;
      }
    }
    memo = malloc(sizeof(int) * (long)nsets * nsets); how = malloc(sizeof(int) * (long)nsets * nsets);
    memset(memo, -1, sizeof(int) * (long)nsets * nsets);
    int c = cost(0, 0);
    printf("F %d d %d s %d bound %d sets %d%s\n", F, d, sF, c, nsets, c < sF ? " VIOLATION" : "");
    if (c < sF) { nviol++; printf("CERT %d\n", F); cert(0, 0); printf("END\n"); }
    free(memo); free(how);
    fflush(stdout);
  }
  printf("DONE tested %d violations %d\n", ntest, nviol);
  return 0;
}
