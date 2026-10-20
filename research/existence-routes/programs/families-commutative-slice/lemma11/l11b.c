/* exhaustive check of Lemma 11 over F_2: no nonsingular n x n matrix has all principal submatrices of order <= s+2
   of rank <= s (n >= s+1).  Matrices as n*n bits; rank by elimination on row bitmasks. */
#include <stdio.h>
#include <stdlib.h>
static int rank_rows(unsigned r[], int k) {
  int rk = 0;
  for (int bit = 0; bit < 8 && rk < k; bit++) {
    int p = -1; for (int i = rk; i < k; i++) if (r[i] >> bit & 1) { p = i; break; }
    if (p < 0) continue;
    unsigned t = r[p]; r[p] = r[rk]; r[rk] = t;
    for (int i = 0; i < k; i++) if (i != rk && (r[i] >> bit & 1)) r[i] ^= r[rk];
    rk++;
  }
  return rk;
}
int main(int argc, char **argv) {
  int n = atoi(argv[1]), s = atoi(argv[2]); int ord = argc > 3 ? atoi(argv[3]) : s + 2;
  long long total = 1LL << (n * n), nonsing = 0, bad = 0;
  for (long long x = 0; x < total; x++) {
    unsigned row[8];
    for (int i = 0; i < n; i++) row[i] = (x >> (i * n)) & ((1u << n) - 1);
    unsigned tmp[8]; for (int i = 0; i < n; i++) tmp[i] = row[i];
    if (rank_rows(tmp, n) < n) continue;
    nonsing++;
    int ok = 1; /* ok = all principal of order <= s+2 have rank <= s */
    for (unsigned T = 1; T < (1u << n) && ok; T++) {
      int k = __builtin_popcount(T); if (k > ord) continue;
      unsigned sub[8]; int c = 0;
      for (int i = 0; i < n; i++) if (T >> i & 1) {
        unsigned v = 0; int cc = 0;
        for (int j = 0; j < n; j++) if (T >> j & 1) { if (row[i] >> j & 1) v |= 1u << cc; cc++; }
        sub[c++] = v;
      }
      if (rank_rows(sub, k) > s) ok = 0;
    }
    if (ok) { bad++; if (bad < 5) printf("COUNTEREXAMPLE n=%d s=%d x=%llx\n", n, s, x); }
  }
  printf("n=%d s=%d over F_2: %lld nonsingular matrices, %lld with all principal submatrices of order <= %d of rank <= %d\n", 
         n, s, nonsing, bad, ord, s);
  return 0;
}
