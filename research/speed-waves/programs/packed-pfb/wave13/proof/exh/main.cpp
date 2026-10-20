// exhaustive: all 2^28 inputs (a, b 8-bit, c 12-bit; c = -2048 is mapped to -2047 inside both, as the ROM never holds it)
#include "Vtop.h"
#include <cstdio>
int main() {
  Vtop *m = new Vtop; unsigned long long n = 0, bad = 0, badx = 0;
  for (unsigned c = 0; c < 4096; c++) for (unsigned a = 0; a < 256; a++) for (unsigned b = 0; b < 256; b++) {
    m->a = a; m->b = b; m->c = c; m->eval(); n++;
    if (m->pa_n != m->pa_r || m->pb_n != m->pb_r) { if (bad < 3) printf("MISMATCH a=%u b=%u c=%u\n", a, b, c); bad++; }
    if (m->pa_x != m->pa_r || m->pb_x != m->pb_r) badx++;
  }
  printf("vectors %llu  new-vs-ref mismatches %llu  broken-vs-ref mismatches %llu\n", n, bad, badx);
  return bad != 0;
}
