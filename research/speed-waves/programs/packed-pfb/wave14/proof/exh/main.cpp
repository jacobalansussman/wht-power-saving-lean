// Exhaustive check: all a, b in [-128, 127] and all 11-bit c in [-1024, 1023]. Reference: ya = a*c, yb = b*c.
// A coefficient is LEGAL for one tap when 127P + 128N <= 2^16 - 1 and 128P + 127N <= 2^16 (i.e. -511 <= c <= 512).
#include "Vtop.h"
#include <cstdio>
#include <cstdint>
static inline long sx(uint32_t v, int w) { return (v >> (w - 1)) & 1 ? (long)v - (1L << w) : (long)v; }
int main() {
  Vtop *m = new Vtop;
  long bad_legal = 0, bad_illegal = 0, n_legal = 0, n_illegal = 0, bad_broken_legal = 0;
  int cmin_bad = 99999, cmax_bad = -99999;
  for (int c = -1024; c < 1024; c++) {
    long P = c > 0 ? c : 0, N = c < 0 ? -c : 0;
    bool legal = (127 * P + 128 * N <= 65535) && (128 * P + 127 * N <= 65536);
    for (int a = -128; a < 128; a++) for (int b = -128; b < 128; b++) {
      m->a = a & 255; m->b = b & 255; m->c = c & 2047; m->eval();
      bool ok = sx(m->ya, 17) == (long)a * c && sx(m->yb, 17) == (long)b * c;
      bool okb = sx(m->ya_br, 17) == (long)a * c && sx(m->yb_br, 17) == (long)b * c;
      if (legal) { n_legal++; if (!ok) bad_legal++; if (!okb) bad_broken_legal++; }
      else { n_illegal++; if (!ok) { bad_illegal++; if (c < 0 && c > cmax_bad) cmax_bad = c; if (c > 0 && c < cmin_bad) cmin_bad = c; } }
    }
  }
  printf("legal coefficients (-511..512): %ld vectors, mismatches %ld %s\n", n_legal, bad_legal, bad_legal ? "FAIL" : "PASS");
  printf("broken copy (correction bit P[15]) on legal coefficients: mismatches %ld %s\n", bad_broken_legal, bad_broken_legal ? "(FAILS as it must)" : "(NOT CAUGHT)");
  printf("illegal coefficients: %ld vectors, mismatches %ld (failing c nearest the legal range: %d and %d) -> the bound is tight\n", n_illegal, bad_illegal, cmax_bad, cmin_bad);
  delete m; return 0;
}
