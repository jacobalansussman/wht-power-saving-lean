// ap0lib2.c -- ap-0: libap0 v2 for Python: 512 / 4,096 through smegi.h (in-kernel interleaving), 1,024 / 8,192
// through the decimation-in-time path (ap0p.c, LD4D, NEON butterfly), 2,048 / 16,384 through the radix-4 conversion
// path with 131,072-point chunks (ap0lib.c).
#include "ap0lib.c"
#include "ap0p.c"
typedef struct { ap0_t *a; ap0p_t *p; } ap0v2_t;
void *ap0v2_new(long N) { ap0v2_t *h = calloc(1, sizeof *h); long m = N; int P = 0; while (m % 8 == 0 && m > 1) { m /= 8; P++; }
  if (m == 2 && P >= 3) h->p = ap0p_new(N); else h->a = ap0_new(N); return h; }
void ap0v2_fft(void *hv, const double *x, double *y, long B) { ap0v2_t *h = hv; if (h->p) ap0p_fft(h->p, x, y, B); else ap0_fft(h->a, x, y, B); }
