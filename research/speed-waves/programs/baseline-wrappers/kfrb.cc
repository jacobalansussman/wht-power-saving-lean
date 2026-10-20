// kfrb.cc -- ar-0: ctypes entry for batched KFR complex<double> rows (numpy layout)
extern "C" { void *kd_new(long n); void kd_exec(void *h, const double *in, double *out);
void kb_batch(void *h, long n, const double *x, double *y, long B) { for (long t = 0; t < B; t++) kd_exec(h, x + 2 * t * n, y + 2 * t * n); } }
