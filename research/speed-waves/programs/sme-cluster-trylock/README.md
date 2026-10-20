# A per-cluster try-lock for the matrix unit in multi-process code (section 3e)

The M4 has one matrix unit per cluster of four P-cores. `mpcore.h` is a multi-process driver that runs W worker
processes calling one kernel in a loop and compares policies for who gets the unit: all NEON; every worker on the
matrix unit; a static split; a try-lock token in shared memory with a NEON fallback (`cauto`, the rule that works);
and a resident owner process. Result: the try-lock gives -17.6 % / -18.9 % on a 4-worker Whisper log-mel and -11 % /
-12 % on a 4-worker double-precision FFT against all-NEON; every worker on the unit is 19 % to 41 % slower; the owner
process fails. Status: checked by computer, one program, not repeated by anyone (read by the fairness critic only).

| file | what it does |
|---|---|
| `mpcore.h` | the driver and the policies (`neon`, `sme`, `static:K`, `token:K`, `cauto`, `owner`) |
| `mp_lm.c` | the workload "Whisper log-mel, 30-second chunks" (kernel in `logmel/`) |
| `mp_fft.c` | the workload "double-precision FFT rows" (kernel from [../sme-fp64-batched-fft](../sme-fp64-batched-fft/)) |
| `cpunum.c` | shows which core a process runs on |
| `logmel/` | the Whisper log-mel kernel of wave 13 (`lm.c`, the generated `kern.h`, `kern_r.h`, `mel_gen80.h`; `gen_mel.py` writes the mel table) |
| `results/` | original outputs (`t8.txt` the final comparison, `t9.txt` against FFTW double, `t10.txt` a repeat) |

```sh
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$FFTW/include mp_lm.c $FFTW/lib/libfftw3f.a -framework Accelerate -o mp_lm
clang -O3 -mcpu=apple-m4 -Wno-deprecated-declarations -I$PFFFT/include -I$FFTW_DOUBLE_NEON/include mp_fft.c \
  kfrd.o <KFR objects> $FFTW_DOUBLE_NEON/lib/libfftw3.a $PFFFT/lib/libpffft.a -framework Accelerate -lc++ -o mp_fft
BREAK=1 R=1 ./mp_lm 1:sme 1:owner          # broken copy: must FAIL
./mp_lm 4:neon 4:cauto 4:sme 4:owner       # one row per policy; RELAXED=1 for relaxed atomics
NEONLIB=kfr ./mp_fft 4:neon 4:cauto
```

Needs: an Apple M4, FFTW single (log-mel), PFFFT, FFTW double with NEON and KFR (FFT). Time: a minute per policy row.
Memory under 1 GB. Check what else runs on the P-cores before timing: one busy core changes every four-core figure.

Not re-run for this revision.
