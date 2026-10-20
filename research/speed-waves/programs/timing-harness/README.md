# Timing harness

`bench.h` is the header-only timing harness that every processor timing of the speed-first waves used.

- Monotonic clock (`CLOCK_UPTIME_RAW` on macOS, `CLOCK_MONOTONIC_RAW` elsewhere).
- For each program it calibrates a call count so that one repeat lasts at least 0.2 s, then runs at least 9 repeats.
- **The programs being compared take turns inside one process**, and the starting program rotates each round, so no
  program is always first.
- Every repeat is printed; best and median are reported in ns per unit (per point by default).
- Two buffer modes: `BENCH_SAME` (one array again and again: warm cache) and `BENCH_SEPARATE` (a ring of arrays of at
  least 64 MiB in total, so each call meets cold data).
- Reference-check helpers to call before timing; on macOS it asks for the user-interactive quality-of-service
  class, which keeps the thread on a performance core.

`bench_example.c` shows the use on a small FFT-like loop:

```sh
clang -O3 -mcpu=apple-m4 bench_example.c -o bench_example
./bench_example 10        # log2 of the size
```

Needs: a C compiler; nothing else. Time: a few seconds.
