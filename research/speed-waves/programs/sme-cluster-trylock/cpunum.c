/* cpunum.c -- does TPIDRRO_EL0 / TPIDR_EL0 carry the CPU number on this macOS, and which numbers are P / E? */
#include <pthread.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
static inline uint64_t ro(void) { uint64_t p; __asm__ volatile("mrs %0, TPIDRRO_EL0" : "=r"(p)); return p; }
static inline uint64_t rw(void) { uint64_t p; __asm__ volatile("mrs %0, TPIDR_EL0" : "=r"(p)); return p; }
static void *th(void *a) {
  int bg = (int)(long)a;
  pthread_set_qos_class_self_np(bg ? QOS_CLASS_BACKGROUND : QOS_CLASS_USER_INTERACTIVE, 0);
  int seen[64] = {0};
  for (long i = 0; i < 200000000; i++) { volatile int x = 0; (void)x; if ((i & 0xfffff) == 0) { uint64_t r = rw(); seen[r & 63]++; } }
  printf("%s: TPIDR_EL0 low bits seen:", bg ? "background" : "user-interactive");
  for (int k = 0; k < 64; k++) if (seen[k]) printf(" %d(%d)", k, seen[k]);
  printf("  | ro=%#llx rw=%#llx\n", (unsigned long long)ro(), (unsigned long long)rw());
  return NULL;
}
int main(void) { pthread_t a, b; pthread_create(&a, NULL, th, (void *)0); pthread_create(&b, NULL, th, (void *)1); pthread_join(a, NULL); pthread_join(b, NULL); return 0; }
