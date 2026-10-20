# Speed-first waves 11-15: real processors and FPGAs (11-13 October 2026)

The rest of this repository counts operations. This section is about time and chip area. On the evening of 11
October 2026 I changed the goal for two days. I wanted algorithms that make real Walsh-Hadamard, FFT and
convolution work **significantly faster on processors, or smaller and faster on FPGAs**, with half of the effort on
each side. The reason was a result I had just published: my circuits with fewer operations were *not* faster inside
a tuned library (see [wht-fewer-operations](https://github.com/jacobalansussman/wht-fewer-operations)). So I went
looking for speed where it actually lives.

I did this work with a team of AI agents (Claude) that I directed. I set the goals and the bar, chose which leads to
push, and decided at each stage what to do next. The agents wrote the kernels and the hardware generators, ran the
measurements, and tried to refute one another's numbers. **Nothing in this section is in Lean.** Every processor
number is a measurement on one Apple M4 Mac mini, except one check that ran on GitHub's Arm machines. Every FPGA
number comes from the open tools (Yosys and nextpnr), not from vendor tools or a board. No person has repeated any of
these measurements. Where I say "a critic", I mean an AI agent that tried to break the claim with its own baseline
and its own programs. It is not a human referee.

**In one paragraph.** Over five waves, about forty developments and two dozen critics' reports, exactly one result
moved a whole application by more than 10 % on every core of the chip. It is a faster number-theoretic transform
(NTT) inside Microsoft SEAL, a homomorphic-encryption library. A loop of real CKKS operations ran 21.2 % faster,
with bit-identical ciphertexts. I then ported it cleanly into SEAL 4.6.2 as an optional build setting that passes
SEAL's whole test suite. On GitHub's Arm server machines the gain is real but smaller. Everything else that survived
is a kernel gain: double-precision FFTs, a Fourier-optics step, Welch spectra and convolution on the M4's matrix unit
(11 % to 87 % less time on one core, depending on the size), and a few FPGA designs that are 20 % to 60 % smaller
with the open tools. These gains shrink, often to nothing, once the whole chip or the whole application is counted. Along the way, a
one-hour probe found that the FFTW most Arm users install has no vector code at all and runs 1.3 to 2.6 times slower
than it could. The section also lists every lead that failed, with its exact scope, so that nobody has to try it
again blind.

## Contents

- [The question and the bar](#the-question-and-the-bar)
- [Words used here](#words-used-here)
- [How the waves worked](#how-the-waves-worked)
- [What I found](#what-i-found)
  1. [Four lessons](#1-four-lessons)
  2. [A faster NTT inside Microsoft SEAL: the one whole-application result](#2-a-faster-ntt-inside-microsoft-seal-the-one-whole-application-result)
  3. [Double-precision Fourier work on the M4's matrix unit](#3-double-precision-fourier-work-on-the-m4s-matrix-unit)
  4. [Single-precision kernels on the matrix unit](#4-single-precision-kernels-on-the-matrix-unit)
  5. [The FFTW that people install on Arm has no vector code](#5-the-fftw-that-people-install-on-arm-has-no-vector-code)
  6. [FPGA results](#6-fpga-results)
  7. [Leads parked because their share of real work is too small](#7-leads-parked-because-their-share-of-real-work-is-too-small)
  8. [What failed, was refuted or fell below the bar](#8-what-failed-was-refuted-or-fell-below-the-bar)
  9. [What was claimed and what survived](#9-what-was-claimed-and-what-survived)
  10. [Traps for anyone measuring on these chips and tools](#10-traps-for-anyone-measuring-on-these-chips-and-tools)
  11. [Wave 15 and the crash](#11-wave-15-and-the-crash)
- [What is still open, with a first experiment for each](#what-is-still-open-with-a-first-experiment-for-each)
- [Programs and data](#programs-and-data)
- [Credits](#credits)
- [Sources](#sources)

## The question and the bar

Which algorithms make real WHT, FFT and convolution uses significantly faster on today's processors, or smaller at
the same speed on FPGAs? I fixed the bar before the first wave:

- **Processors:** at least 10 % less time than the strongest available baseline, like for like, at sizes real users
  run.
- **FPGAs:** at least 10 % less area at equal throughput and equal clock.
- **Or** a capability that no baseline has, with a measured advantage.

"Strongest available baseline" carries most of the weight. Most claims in this section first beat a weak baseline
and then shrank when a critic built a better one. Section 9 lists every such correction.

## Words used here

- **M4**: the Apple M4 Mac mini that every processor timing comes from (4 performance cores, 6 efficiency cores).
  **M1**: a second Mac, which has no matrix unit; it was asleep and unreachable for most of these waves.
- **P-core / P-cluster**: a performance core / the four of them together.
- **SME**: the M4's matrix unit (Arm's Scalable Matrix Extension). There is **one per cluster**, shared by the four
  P-cores, so four threads cannot each have their own. It is reached through a special "streaming" mode.
- **NEON**: the ordinary 128-bit vector instructions that every 64-bit Arm core has.
- **fp32 / fp64**: single / double precision. **complex128**: double-precision complex numbers, numpy's default.
- **In cache / from DRAM**: the data is already in the processor's cache / comes from main memory.
- **Baseline libraries**: FFTW, KFR, PFFFT and pocketfft (open FFT libraries; pocketfft is the one inside numpy and
  scipy), Apple's vDSP (part of Accelerate), and FFHT (a fast Walsh-Hadamard library).
- **NTT**: the number-theoretic transform, an FFT over integers modulo a prime. It dominates homomorphic encryption.
  **CKKS** and **BFV** are the two encryption schemes timed here.
- **P&R**: placed and routed with the open FPGA tools. **Synthesis only**: synthesised but never placed and routed,
  so there is no clock and no placed area. FPGA clocks are the tools' estimates, never a board measurement.
- **LUT4, slices, cells**: FPGA area units. **DSP**: an FPGA's hard multiplier block.
- **iCE40-UP5K**: a small, battery-class Lattice FPGA. **ECP5, CrossLink-NX**: mid-range Lattice FPGAs.
  **Xilinx 7-series, UltraScale**: AMD's FPGA families (only synthesis here).
- **PFB**: polyphase filter bank, the channeliser at the front of a radio telescope.
- **SATD**: the Hadamard-based cost that video encoders compute millions of times per frame.
- **Fair rival**: the obvious design written well, without the developer's new trick. From wave 14 on, the critics
  built one for every FPGA claim.

## How the waves worked

Each wave had eight developments, four on processors and four on FPGAs, run in parallel for about 80 minutes each.
A recombination step then tabulated the results, formed hybrids and picked the top claims. One critic per top claim
tried to refute it with a stronger baseline; a fairness critic read the whole wave; two "pool hunters" collected new
ideas. A note-writer then kept a running list of live leads, parked leads, closed doors and corrections. Waves 11 to
14 ran overnight from 11 October 21:35 to 12 October 06:40 (US Eastern time). Wave 15 was cut short by a machine crash
(section 11). Four scouts measured, between waves 11 and 12, **where the transform is most of the work** in real
software and hardware, because wave 11 had shown that kernel gains on small shares are worthless.

Measurement rules, applied to every number below:
- **Timing**: one shared harness ([programs/timing-harness](programs/timing-harness/)). At least 9 repeats of at
  least 0.2 s each; the programs being compared take turns inside one process, so none is always first; every
  repeat is printed; best and median are reported. Timings ran on a P-core under a machine-wide timing lock. The
  other agents kept the efficiency cores busy, with load averages of 10 to 90, which is why critics re-ran the
  important figures.
- **Correctness before timing**: every kernel was checked against a double-precision or exact reference first, and
  a deliberately broken copy had to fail the same check.
- **FPGAs**: one measuring tool for all chips ([programs/fpga-measure](programs/fpga-measure/)). The design sits in
  a small harness so that its ports do not need pins; 5 placement seeds; the median is reported. Bit-exact
  simulation against a bit-true model, and a broken copy that must fail.

## What I found

Each finding carries its status in words. For measurements I use them like this:
- **checked by computer**: a program measured it, and an independent AI critic re-measured it with its own baseline
  and programs. I say how many critics.
- **checked by computer, one program, not repeated by anyone**: measured once, and no critic re-ran it.
- **refuted**: a later measurement broke it. I say by what.
- **open**: not settled. Each open item has a first experiment at the end.

Where notes disagree, the later critic's figure is the one given, and I say so.

### 1. Four lessons

These are the most reusable results of the section. Each one was learned by losing a claim.

1. **Kernel gains are rarely application gains.** Waves 11 and 12 had no result that moved a whole application by
   10 %. The SEAL NTT (section 2) is still the only one: it moves a real loop by 21 % because the NTT is 59 % to 66 %
   of that loop. Measure the transform's share of the whole job first.
2. **Every headline shrank under a stronger baseline.** The stronger baseline was a better library (FFTW double
   precision built with NEON *and* a working timer, or KFR), more cores (one matrix unit per cluster: on all four
   P-cores most matrix-unit gains vanish), colder caches, or a different design.
3. **Credit the known part.** Most FPGA wins were a known mechanism in a new use: Xilinx's WP486 packing, vertical
   common-subexpression sharing, LUT-RAM transposition, carry-cell absorption. Say which part is new and measure that
   part alone.
4. **The obvious rewrite takes a large share of the gain.** Plainly written gated RTL took a fifth of one FPGA FFT's
   gap; DSP code in the vendor-manual style took two-thirds of another's; letting the tool map a plain `*` to a DSP
   block cut one claim from -82 % to -56 %; a 6-line batched-FFT Welch in Python cut another from -85 % to -62 % /
   -66 %. Build the fair rival before claiming.

### 2. A faster NTT inside Microsoft SEAL: the one whole-application result

**The idea.** SEAL's NTT on Arm is scalar 64-bit integer code. A scout measured that the NTT is 67 % to 68 % of
SEAL's CKKS work on the M4. Primes below 2^50 allow exact NTT arithmetic in double-precision floating point with
fused multiply-add: every value stays a whole number small enough to be held exactly in a double, and the rounding of
each product is recovered exactly. This is not my invention. J. van der Hoeven and G. Lecerf (2024) and FLINT's
`fft_small` already use this arithmetic. What is new here is putting it inside SEAL, measuring what it does to real
SEAL operations, and making it a clean, tested build option.

**Wave 13: the first measurement.** A NEON kernel of this kind replaced SEAL 4.1.2's NTT for primes below 2^50.
- One NTT of 2^14 points with a 50-bit prime: 39.6 against 61.6 µs (-35.8 %).
- A loop of real CKKS operations at N = 2^14 (moduli of 60, 7 × 40 and 60 bits): **93.6 against 118.8 ms
  (-21.2 %)**; -14.4 % against a tuned scalar control; **the same -21.2 % with four loops on the four P-cores**,
  -22 % with eight. The ciphertexts were word-for-word equal to SEAL's, and SEAL's 276 tests passed. In this loop the two 60-bit primes
  ran a tuned scalar radix-4 NTT, also written in wave 13; with SEAL's own code for them the gain was -19.1 % (the
  developer's run, not repeated by the critic).
- SEAL's own BFV default moduli at N = 2^12 to 2^14 are all below 2^50: **-16.3 % to -20.8 %**. The developer had
  first said "default moduli only -8 %"; the critic found that wrong in the good direction.
- CKKS at N = 2^15 (60, 10 × 40, 60): -27.0 %. The tutorial size N = 2^13 (60, 40, 40, 60): -10.8 % (-7.1 % against
  the scalar control), because half of its NTTs are on 60-bit primes, which keep the old code.

Status: **checked by computer**, one critic re-measured it and a fairness critic read it.

**Wave 14: primes of 55 to 60 bits.** Two routes were tried for the large primes that the floating-point method
cannot take. A radix-4 scalar NTT with one overflow guard per pass reached -16.7 % on BFV at N = 2^15 against SEAL as
shipped, but only **-6.3 % against our own scalar control: below the bar**. An exact integer route on the matrix unit
(16-bit limb products) was **+29 % to +49 % slower**: moving partial sums out of the matrix unit and the final
87-bit reduction cost more than the products saved. Status: **checked by computer, one program, not repeated by
anyone** (read by the fairness critic only). Both are closed doors.

**12-13 October: a clean port into SEAL 4.6.2.** I had the kernel rebuilt as an optional part of SEAL 4.6.2:
- One new build setting, `SEAL_USE_NEON_NTT`, **off by default**. When on, every prime below 2^50 uses the NEON code;
  larger primes and transforms of fewer than four coefficients keep SEAL's code. The code plugs into the same four
  functions where SEAL already plugs in Intel's HEXL library for x86.
- It only takes effect on 64-bit Arm with GCC or Clang; elsewhere CMake switches it off with a message. With it
  off, SEAL's compiled library objects are byte-for-byte those of unmodified SEAL 4.6.2.
- The new file always compiles with `-fno-fast-math -ffp-contract=off` and refuses to compile under fast-math. At
  run time it is used only under round-to-nearest; the transforms leave the floating-point status flags as they found
  them (building the tables does not; see the skeptics' items below); it uses only baseline 64-bit Arm instructions (no Apple-specific code, no matrix unit). A comment at the top of the
  file gives the exactness proof with its bounds. The proof is argued by hand and is not checked in Lean.
- Fully reduced results are identical to SEAL's. SEAL's "lazy" variants may return any value below 4p or 2p; the new
  code returns the fully reduced value there, which is congruent and in range (HEXL is allowed to do the same) but
  not always the same 64-bit word.
- Each NTT table grows by 16 bytes per coefficient. Start-up memory of SEAL's benchmark program went from 467 to
  502 MB (+7.5 %).

Tests and checks:
- **SEAL's full test suite passes**: 372 of 372 tests with the setting off (SEAL's 368 plus 4 new ones) and 373 of
  373 with it on, in Release and Debug builds on the M4, and with GCC 13.3 and Clang 18.1 on GitHub's Arm machines.
  On an x86 machine the setting switched itself off and all 372 tests passed.
- **The new tests catch real mistakes**: of 12 deliberately planted bugs, they catch 11. The survivor tells the
  reduction planner that rounding errors are 6 times smaller than they can be; random inputs almost never reach that
  worst case, so for that margin the written proof protects the code, not the tests.
- **Two independent AI skeptics** then tried to break the port. **Neither refuted it.** They compared about 36.6
  million transforms against SEAL's scalar code with 0 mismatches; found CKKS, BFV and BGV results byte-identical
  with the setting on and off in all four rounding modes; and checked the reduction planner with exact rationals in
  7,130 cases. They left small items open: a floating-point "inexact" trap, if a program enables it, would stop the
  process; building the tables raises the inexact flag; GCC needs two more flags to rule out unsafe math per flag;
  an empty header is installed when the setting is off; and SEAL's own continuous integration tests only x86, so
  SEAL would need an Arm job.

Status of the port: **checked by computer** (SEAL's tests on two chips and two compilers, planted bugs, two
skeptics).

**Speed of the port, SEAL's own benchmark program, M4** (generic 64-bit Arm build; 7 off/on pairs; SEAL's 128-bit
default primes; n = 4,096 to 16,384):

| operation | change with the setting on |
|---|---|
| the NTT itself | -35 % to -58 % |
| relinearisation and rotation (key switching) | -24 % to -30 % |
| CKKS rescale | -30 % to -33 % |
| BFV decryption | -20 % to -26 % |
| BFV multiplication and squaring | -5 % to -9 % |
| multiply then relinearise (the sum of the two measured operations) | BFV -7 % to -12 %, CKKS -21 % to -26 % |
| CKKS multiply and add (no NTT inside; controls) | no change |
| everything at n = 32,768 (primes of 55-56 bits; the old code runs) | no change, except BFV batch encoding and decoding (-44 %, -53 %), which run on a 20-bit plaintext prime |

A control experiment ran the same program with the NEON code switched off at run time, so that code layout and
memory were identical. It matched the table on the rows with an NTT inside (a median gap of 0.8 points, at most 4.6). It also
showed that an apparent -20 % on CKKS decryption, which has no NTT inside, is not the port: most likely it is where
the machine code happened to sit in memory (an inference from the control runs; the loop itself was not traced).
The BFV decryption figure above is the control's. Status: **checked by computer, one program** (the skeptics checked
correctness, not these timings).

**Speed on GitHub's Arm machines** (Arm Neoverse N2, 2 virtual processors, two full runs on fresh machines that
agree within about half a point; 7 pairs each):
- The fully reduced forward NTT: **-27 % to -30 %** with both compilers (-55 % to -57 % on the M4).
- Public-key encryption: -8 % to -15 %; CKKS public-key encryption clears the bar with both compilers.
- Relinearisation and rotation: **about -4 % to -5 % with GCC and about -9 % to -10 % with Clang** (-23 % to -28 % on the M4).
  These follow the lazy forward NTT, which SEAL's own code already runs almost as fast as the new code on this core.
- BFV multiplication: -1 % to -3 %.
- With GCC at n = 32,768, SEAL's own fully reduced forward NTT ran 3 % to 4 % slower in the build with the setting
  on. The machine code of that loop is identical; only its alignment moved. Clang shows nothing of the kind.

Status: **checked by computer, one program** (two runs, two compilers). Why the two chips differ was not measured.
A plausible reason is the number of vector units per core, but that is a guess, not a measurement.

**What this means.** On the M4 the port is a large, real gain for every SEAL user whose primes are below 2^50,
which includes SEAL's own defaults at the three middle sizes. On the one non-Apple Arm design I could test, it is a
clear gain for the fully reduced transforms and encryption and a small one for key switching. **I am holding the
patch itself back** from this revision until I decide whether and how to offer it to SEAL's maintainers; the result
tables are in [data/seal-ntt-port](data/seal-ntt-port/). Not done: a run on the M1; a build with Intel's HEXL; any
chip other than the M4 and the Neoverse N2; a port to OpenFHE (begun in wave 15, never finished).

### 3. Double-precision Fourier work on the M4's matrix unit

A wave-11 probe found that the matrix unit does double-precision outer products at about 248 billion
multiply-adds per second, 8.3 times NEON's double-precision rate, and that no library used it. A scout found that
double precision, numpy's default, is the worst-served precision on the Mac, and that in Fourier optics and in
scipy's convolution routines the FFT is 91 % to 97 % of the work. Five leads followed.

**3a. Batched complex128 FFT.** A transform of 8^P points as P radix-8 passes (three for 512 points, four for
4,096), each pass a complex 8 × 8 DFT matrix applied with outer products; 1,024 and 2,048 points add one NEON pass;
numpy's memory layout. One P-core, batch 64:
- In cache: **-41 % at 512 points** (0.688 against 1.171 ns per point for KFR, the strongest library there) and
  **-44 % at 4,096 points** (0.952 against 1.709 for FFTW double with NEON). The 4,096 figure was first -52 %; a
  later rebuild of the strongest baselines found that the FFTW plan used before was 17 % to 20 % slow, and the later
  figure is the one given.
- From DRAM: -50 % at 512 and -52 % at 4,096.
- At 1,024 and 2,048 points (one NEON pass plus the matrix unit) only -15 % (1,024, in cache) and -16 % (2,048, from
  DRAM).
- On the whole P-cluster (one matrix-unit thread plus three NEON threads against four KFR threads): -26 % at 512
  and -18 % at 4,096.
- From Python through ctypes: -23 % and -26 % against `scipy.fft` with 4 workers.
- Error 1.5 to 1.8 times the libraries' (relative RMS 4.5e-16 to 5.5e-16).

Status: **checked by computer**, one critic and the baseline rebuild. Its correctness check was re-run for this
revision from the published copy and gave the same errors.

**3b. A Fourier-optics propagation step** (`ifft2(fft2(E) * H)` in complex128, as in angular-spectrum and
beam-propagation codes). One P-core, grids from DRAM, against FFTW double with NEON and a working timer using its own
2-D plan (in place, PATIENT) or KFR, whichever was faster:
- **-38 % / -47 % / -62 %** at 1,024² / 2,048² / 4,096² (5.0 against 8.1 ms at 1,024²), and -49 % at 512².
- The output is in index-reversed order. Natural order costs extra; with it the gains are about -30 % / -42 % / -56 %
  (worked out from measured parts, not timed as one).
- **On the whole P-cluster the 1,024² claim is refuted: it is a tie** with FFTW's own 4-thread 2-D plan. At 2,048²
  and 4,096² the cluster gain is -12.7 % and -20.7 % against four NEON threads running the same strip plan.
- Error 1.2e-15 to 2.1e-15. If single precision is acceptable, the M4's GPU (through MLX, complex64) is 2.6 to 3.7
  times faster than one matrix-unit core. This result is for users who need complex128.

Status: **checked by computer**, one critic and the fairness critic; the cluster claim at 1,024² **refuted** by the
same critic.

**3c. A one-pass double-precision Welch spectrum.** Window, FFT and |X|² accumulation fused in one pass. One P-core,
2^22 samples, 50 % Hann overlap, the whole call, against FFTW double with NEON and the real-input routes:
- **-33 % / -27 % / -11 % / -41 %** at segment lengths 512 / 1,024 / 2,048 / 4,096 (4.59 against 6.87 ms at 512;
  5.65 against 9.55 ms at 4,096). The developer's -30 % and -15 % at 1,024 and 2,048 shrank under the critic.
- **On four P-cores no size clears the bar** (-1 % to -10 %): the developer's -10.3 % and -11.2 % are **refuted**;
  the 1,024 claim had only vDSP in the race.
- From Python: -62 % / -66 % against a 6-line batched `rfft` Welch with 4 workers. Most of that is "Welch in C at
  all", which would help every processor.
- 239 checks against scipy, largest error 4.7e-15.

Status: **checked by computer**, one critic and the fairness critic.

**3d. numpy's convolve and correlate on the matrix unit** (wave 15). numpy always computes the direct sum, at 5.6 to
11 billion multiply-adds per second. A Toeplitz-tile kernel in double precision, wrapped as a drop-in for
`np.convolve` and `np.correlate` with numpy's modes and argument rules:
- Whole Python call, double precision, mode 'same', 10^6 samples, against the best of numpy, scipy (`convolve`
  auto, `oaconvolve`, `fftconvolve`) and vDSP's `convD` at each size: **-74 % / -85 % / -87 % / -81 % / -68 % /
  -70 %** at 9 / 33 / 129 / 257 / 513 / 1,025 taps. 10^4, 10^5 and 10^7 samples look the same (-50 % to -86 %
  from 9 to 1,025 taps); complex128 -56 % to -88 % from 3 to 257 taps.
- The kernel runs at 240 billion multiply-adds per second at 1,025 taps, 97 % of the unit's double-precision peak.
- Accuracy: a blocked variant keeps the error within 1.4 times numpy's median (worst seen 3.6 times, all below
  3e-16 of the sum of |terms|) for 0 % to 7 % more time; NaN and infinity inputs give exactly numpy's NaN pattern.
- **Against the strongest FFT route I could build (a double-precision PFFFT overlap-save in C), the direct kernel
  wins only below about 450 taps**: -33 % at 257 taps, +13 % at 513, 1.9 times slower at 1,025. The big numbers
  above are against what numpy and scipy ship, not against the best possible FFT method.
- It does not clear the bar for 3-tap kernels on 10^7 samples (memory-bound), for calls under about 1,000 outputs
  (a fixed cost of about 1.7 µs per call), or for complex128 at 1,025 taps.
- It has no per-cluster guard yet (3e), so several Python processes calling it at once would contend for the one
  matrix unit.

Status: **checked by computer, one program, not repeated by anyone**: wave 15 ended before a critic could read it.
Its correctness checks were re-run for this revision from the published copy: 534 of 534 pass, the broken copy fails
4 of 4, and the NaN pattern is identical.

**3e. A rule for the matrix unit in multi-process code.** Because one matrix unit serves a cluster, several worker
processes that all enter it slow each other down: **19 % to 41 % slower than all-NEON, with the slowest 1 % of calls
1.6 to 2.8 times slower**. A per-cluster try-lock in front of the matrix unit, with a NEON fallback for whoever does
not get it, gives **-17.6 % / -18.9 %** on a 4-worker Whisper log-mel front end (-15 % with 8 workers) and -11 % /
-12 % on a 4-worker double-precision FFT, against all-NEON. A resident "owner" process that serves the others through
shared memory fails: the hand-off alone costs 9 to 25 µs per call. Status: **checked by computer, one program, not
repeated by anyone** (read by the fairness critic only). It is a design rule, not a speed-up on its own.

### 4. Single-precision kernels on the matrix unit

These came from wave 11, the first wave, before the scouts.

- **A Toeplitz-tile FIR filter for long single-channel signals.** fp32, one thread, 2^20 samples, in cache, at least
  about 1,000 outputs per call: **-59 % / -53 % / -50 %** at 16 / 64 / 256 taps against Accelerate's own
  matrix-unit `sgemm` on block-Toeplitz panels. The developer claimed -75 % / -87 % / -68 % against vDSP's `conv` and
  PFFFT; two critics independently built the stronger `sgemm` baseline (a 40-line program, 1.6 to 3.5 times faster
  than the developer's best). From DRAM -42 % / -49 % / -52 %. It **loses badly on short calls** (+127 % to +481 % at
  128 to 256 outputs per call) and above about 500 to 1,000 taps an FFT convolution wins (the critic's estimate from measured parts). Status: **checked by computer**, two
  critics; correctness re-run for this revision with identical output.
- **A batched 256-point complex FFT.** -28 % per core against vDSP's DFT at batch 64; -31 % to -33 % from DRAM; it
  loses at batch 1 (+6.6 %). **On the whole chip only -5 % to -12 % in cache, and nothing from DRAM**: vDSP itself
  uses the same matrix unit. Status:
  **checked by computer**, one critic; correctness re-run for this revision with identical output.
- The Hadamard rotation used by quantised language models and the Whisper spectrogram are in section 7: both are
  fast kernels with a tiny share of their applications.

The NEON Walsh-Hadamard transform written as the shared baseline in wave 11 is itself 8 % to 24 % faster than
FFHT's NEON code at each of the six sizes tested from 2^8 to 2^22 (fp32, one thread). Status: **checked by computer, one program, not
repeated by anyone**; its correctness check was re-run for this revision with identical output.

### 5. The FFTW that people install on Arm has no vector code

A one-hour probe in wave 14 looked at the FFTW builds that most people actually install. conda-forge's FFTW 3.3.10
for `osx-arm64` and `linux-aarch64` contains **0 NEON codelets** in either precision; Homebrew passes its SIMD flags
on Intel only; Debian's arm64 packages have NEON in single precision only. FFTW's own NEON code exists (a correct
build has 583 NEON codelets in each precision). With the same API calls and plans, one P-core, in cache, conda-forge's
library takes **1.63 to 1.68 times as long for complex double-precision transforms, 1.27 to 1.38 times for real
double, and 2.0 to 2.6 times for single precision**.

This is known to be wanted and is still not shipped. conda-forge's `fftw-feedstock` issue #94 disabled NEON after a
benchmark whose reporter later withdrew it: FFTW's `configure` turns on double-precision NEON only when the host CPU
is called `aarch64`, and Apple reports `arm64` (FFTW issue #251). What the probe adds is the size of the loss on
single transforms and a reproducible recipe (section 10).

The same probe found that pocketfft's batched vector path (the FFT inside numpy and scipy) is **30 % slower than its
scalar path** at 64 rows of 4,096 complex128 points (3.86 against 2.96 ns per point, one run). My guess, not measured, is
that it copies two rows into a buffer the size of a P-core's whole L1 cache.

Status: **checked by computer, one program, not repeated by anyone** (no critic read the probe).

### 6. FPGA results

All FPGA numbers come from the open tools. None was checked with Lattice's or AMD's own tools, none on a board, and
there is no power figure anywhere. That is the weakest side of the whole section.

**6a. A telescope filter bank at half the multipliers.** In the SKA-Low station FPGA (a Xilinx KU040), the
channeliser uses 77 % of the used DSP blocks, the chip's binding resource (from the published firmware paper). Two
8-bit signals can share one 27 × 18 DSP48E2 multiply when they share a coefficient. **The halving itself is Xilinx's
known WP486 packing.** What is new is the precision: an exact guard derived from each branch's (or each tap group's)
coefficient sum, an offset through the DSP's C port, and tap groups graded by coefficient size.
- Wave 13: 224 → 112 DSPs per slice of the design. The new part saves about 17,000 LUTs per chip and gains 1.6 to
  2.1 dB of stop-band over the obvious WP486 design. Status: **checked by computer**, one critic, **synthesis only**.
- Wave 14, graded tap groups: **7 instead of 14 DSPs per signal lane with every filter requirement met** (ripple
  ±0.19 dB, alias edge -62.1 dB, far stop-band -83.2 dB, against -81.4 dB for per-tap WP486 packing at ten times the
  extra LUTs). Worked out for the whole chip: -33.4 % of the channeliser's DSPs, -25.6 % of the chip's. On 7-series
  DSP48E1 blocks one variant meets everything at -43 % of the DSPs. Status: **checked by computer, one program, not
  repeated by anyone** (the fairness critic read it), **synthesis only**.
- Unpacking is bit-exact: all 2^28 inputs of one DSP were checked exhaustively, and a broken copy fails on 41
  million of them. This exhaustive check was re-run for this revision with identical output.

**6b. A streaming FFT with its multiplexers folded into carry cells.** A radix-2² single-path delay-feedback FFT, one
sample per clock, with each stage's multiplexers absorbed into the adders' Lattice carry cells (CCU2):
- **About -20 % slices against the best generic design** (the same gated arithmetic written as plain RTL) on
  CrossLink-NX at 256 points (worked out from measured parts), -21 % LUT positions measured; on ECP5 about 4 points
  less. NX clock +10 %. A minimal range-Doppler radar chain: -14.3 % slices, one seed.
- The developer's -52 % against a published open FFT was dropped: unequal clock and block RAM.
- **Only against Yosys.** Lattice patent US 9,390,210 covers absorbing an adder's output multiplexer into the
  adder's LUTs, so Lattice's own tools may already do this. That was not tested.

Status: **checked by computer**, one critic, open tools only. The 20,000 random checks of the four fused carry-cell
functions were re-run for this revision: 0 failures.

Wave 15 began packing two real radar channels into one complex FFT in that chain. The packed and unpacked chains
were built and simulated bit-exact (114,688 outputs, 0 mismatches; a broken split fails), but the area table was never
written. Status: **open**.

**6c. Spectrum engines on the smallest DSP-equipped Lattice part (iCE40-UP5K).**
- **A folded 4,096-point FFT** with all arithmetic in 6 of the chip's 8 DSP blocks: **235 logic cells** plus 6 DSPs,
  8 block RAMs and 2 large RAMs, reproduced on 10 seeds. **-40 % against the best DSP-multiplier rival** (393 cells),
  -56 % against the obvious design (529). The clock is **50 MHz**, the datasheet's cap for a multiply that skips its
  pipeline register (the tool said about 61 MHz); at 48 MHz it does 1.85 million samples per second. Bit-exact in RTL
  (110,592 samples) and in the netlist as far as the critic's simulation got (53,248 samples, 13 of 25 frames).
  Its 16-bit datapath reads a -70 dBFS noise floor about 16 dB too high (21 dB at 16,384 points; a bit-true model,
  one program), and a Cortex-M4F microcontroller could probably do the vibration-analysis job it was aimed at
  (recalled from published figures, not measured). Status: **checked by computer**, one critic. Its precision model was re-run for this revision with identical output.
- **A 25-bit folded 8,192-point FFT** (wave 15) fixes the precision: 274 cells, 6 DSPs, 16 block RAMs, 2 large RAMs;
  the noise floor is read within 0.73 dB (the 16-bit design reads it 33 dB high), for 10 cells more than a 16-bit
  design of the same size. Rivals at the same widths: logic-only 709 cells (-61.4 %), a behavioural
  multiply-accumulate 645 (-57.5 %), twiddles in three DSPs 408 (-32.8 %). Clock about 50 MHz (same cap). Bit-exact
  in RTL on 106,496 samples; the netlist simulation was not finished. Status: **checked by computer, one program, not
  repeated by anyone**.
- **A streaming 4,096-point FFT** that fits the chip at all: 4,849 cells (92 %), 30.9 MHz over 12 seeds, 47.7 dB
  signal-to-noise. It is a capability; its only rivals are ablations of itself. Status: **checked by computer**, one
  critic.
- **A streaming 4,096-point Walsh-Hadamard transform** with adders and multiplexers hidden in the idle DSP blocks:
  253 against 633 cells (-60 %), about 58 MHz. For the slow uses it was aimed at, a folded design is smaller still (189
  cells, no DSP). Status: **checked by computer**, one critic.
- A free improvement for any halving butterfly: feed the bit that the halving drops into the sum's carry-in. It gains
  6 dB of signal-to-noise for 6 to 13 cells; on the UP5K the DSP's carry-in pin does it for free. Status: **checked
  by computer** (bit-true models), one critic.

**6d. Probes, not criticised.** An NTT datapath for homomorphic encryption on FPGA, using the same one-guard radix-4
schedule as section 2: -18.5 % LUT4 and -14.7 % slices for the fabric of one radix-4 group on ECP5 (the DSP products
excluded). A fast-FIR structure across frames in a telescope filter bank: 21 instead of 28 multiplies per two outputs
(-25 %). Status of both: **checked by computer, one program, not repeated by anyone**; their programs are not in this
revision.

**6e. HEVC motion estimation: a capability lead that turned out to rest on a wrong design.** An open HEVC encoder's
fractional motion estimation needs 397 multipliers, and an ECP5-85F has 156. With hand-written shift-add
interpolation filters it placed at 73 % of the chip with 11 multipliers (wave 14), but routing never finished. On
Xilinx (synthesis only), the interpolator tricks saved 18.6 % of the LUTs and 18.2 % of the DSPs in the module against a fair
rival, but only **-7.6 % LUTs on the whole motion estimation: below the bar** (DSPs -17.8 %). Then wave 15 found that
the faster SATD engine used in those figures **is not a drop-in replacement inside the encoder**: with it the motion
estimation writes different motion vectors for 88 % of entries, because its result arrives later than the encoder's
schedule expects, so each cost uses the previous block's SATD. So **the 73 % placement and the Xilinx figures for the whole motion estimation were for an
incorrect design**; the module-level interpolator figures do not depend on the SATD engine. The correct design needs about 83 % of the chip, and 1080p at 30 frames per second needs 67.8 MHz against a pre-route estimate of
30.6 MHz. Status: the finding that breaks the earlier figures is **checked by computer, one program, not repeated by
anyone**; the lead does not clear the bar. Its programs are not in this revision: they live inside the encoder's
source tree, whose licence forbids modification and redistribution.

### 7. Leads parked because their share of real work is too small

Each is alive as a kernel, with the number that survived, but none moves a real job enough to push now. Their
programs are not in this revision.

| lead | surviving number | status | why parked |
|---|---|---|---|
| Hadamard rotation for quantised language models, matrix unit, fp32 | -73 % at 4,096 points, -64 % at 14,336 (28 × 512) against our NEON WHT (at 14,336 with an `sgemm` stage for the factor 28), one thread, in cache; 15.5 times faster than ExecuTorch's shipping 28·2^k code | checked by computer, one critic | about 0.2 % of a decode token; a whole prefill layer chain -10.5 % on one core, a tie with two cores |
| Whisper log-mel front end, matrix unit, fp32 | spectrogram 80.0 against 163.1 ns per frame for FFTW (-50 %); whole log-mel -40 % on one core, -19.5 % to -22.9 % on four | checked by computer, two critics | from audio files only -1.6 % to -3 % (decoding dominates); under 1-2 % of Whisper inference |
| bank of 8 to 64 short FIR filters, matrix unit | -29 % (16 × 512 taps), -31 % (8 × 1,024) at 128-sample callbacks, warm caches only | checked by computer, one critic | loses 29 % to 50 % with cold caches; about 0.1 % of a core |
| streaming SATD with one's-complement abs-max (FPGA) | 4 × 4: -22.9 % (iCE40), -15.9 % (ECP5), -20.5 % (Gowin); 8 × 8 against LUT-RAM transposition: ECP5 -24 % LUT4 but -2 % to +1 % slices | checked by computer, one critic per size | inside the encoder's motion estimation +0.9 % |
| fused last level in Hadamard-code decoders for 5G control channels (FPGA) | (32,11): -31.3 % LUT4 on ECP5 against the strongest unfused decoder | checked by computer, one critic; the critic checked the fused back end's per-clock logic for all inputs with Yosys's `sat` solver (a transcription of the logic; the pipelined Verilog is covered by simulation) | about 0.15 % of a base-station FPGA |
| the same decoder on NEON for software base stations | -18 % for the new part (the -61 % first claimed is any fast transform against a brute-force correlator) | checked by computer, one program (fairness critic read it) | not yet timed against the real base-station code |
| fp16 lanes for rotate-and-quantise | -20.2 % against fp32 at 4,096 points, equal quality | checked by computer, one critic | 0.2-0.5 % of a decode token |
| polar-decoder element in butterfly form (FPGA) | ECP5 -32 % against a plainly written element (not the strongest) | checked by computer, one program (a pool hunter's probe) | decoders are dominated by memory and sorting |
| ternary-adder Hadamard on Xilinx 7 | -13.9 % by Yosys's count, synthesis only | open: the fairness critic called it unproven (counting rule disputed) | needs vendor place-and-route |
| fused last level in a locality-sensitive hash | -10.4 % to -14.6 % against 3 NEON WHTs plus argmax | open: the baseline does not fuse its argmax | needs a stronger baseline |

### 8. What failed, was refuted or fell below the bar

Exact scope matters: each line closes only what it says. Status: **refuted** where marked; otherwise **checked by
computer** (failed or below the bar), by a critic or by the developer; the Sources section maps each line to its report.

| lead | result |
|---|---|
| int16 arithmetic for a Hadamard rotate-and-quantise on NEON | **refuted**: an fp16 kernel with the same plan is as fast or faster at equal quality; the gain was 16-bit lanes, not integers |
| matrix-unit Welch on four P-cores | **refuted**: no size clears (-1 % to -10 %) |
| matrix-unit optics step on the whole P-cluster at 1,024² | **refuted**: a tie with FFTW's own 4-thread 2-D plan |
| exact 55-60-bit NTT through 16-bit limb products on the matrix unit | +29 % to +49 % slower than a tuned scalar NTT; the whole NTT +37 % |
| radix-4 double-precision NTT on NEON | 3 % to 9 % slower than radix 2 at 2^13 to 2^16 |
| one-guard radix-4 scalar NTT for 55-60-bit primes in SEAL | -6.3 % against our own scalar control: below the bar |
| a resident matrix-unit owner process fed through shared memory | 9-25 µs hand-off and 22-25 µs copy per call; 0 % to 8 % slower on FFTs |
| every worker entering the matrix unit without a guard | 19 % to 41 % slower than all-NEON |
| smaller streaming chunks for the four-core matrix-unit Welch | every mix slower: the shared unit is the limit, not cache size |
| double-precision FFTs or Welch through Accelerate's BLAS (`zgemm` per pass) | 4 to 7 times slower than NEON libraries |
| a Python four-step for one long complex128 FFT from numpy's batched calls | no fix (0.95 to 1.57 times numpy's single call) |
| matrix-unit FIR above about 1,000 taps; short calls of 128-256 outputs; a Karatsuba layer; a second matrix-unit thread with warm data | FFT 2.7 to 3.2 times faster at 4,096 taps; +127 % to +481 %; slower; 0 % |
| matrix-unit FFT at batch 1; unpipelined 320-512 points in cache; twiddles folded across the batch | +6.6 %; tie or loss; slower |
| bf16 or int16 outer products on the M4; double precision emulated from int8 | no faster than fp32; 71 against 248 billion multiply-adds per second |
| Hadamard as Accelerate `sgemm` (H64·X·H64) | ties NEON: the 3-pass plan is the gain, not the matrix unit |
| matrix-unit rotation of a single 256-point vector; fp16 input and output through NEON conversions | +12 %; +3 % to +52 % |
| matrix-unit rotation chained with an int8 matrix multiply in a language-model layer | -10.5 % on one core, a tie with a second core: the rotation's output cannot stay in the matrix unit |
| a matrix-unit direct head in one long partitioned convolution | +8 % to +19 % |
| `sgemm` block-Toeplitz FIR at 64-128-sample callbacks | 3 to 60 times slower than an FFT convolver |
| more matrix-unit work in the Whisper log-mel (mel product inside the unit) | -5 % to -8 % on one core, -1.3 % to -2.5 % on the cluster: below |
| binaural rendering on the matrix unit | -9 % to -11 % on one core at 128-sample callbacks; slower on four cores and from DRAM |
| SATD with int8 matrix or dot-product instructions in x264 | 12 % to 13 % slower than x264's NEON assembly one-to-one; a shared-source form gains only 0.2 % to 0.4 % frames per second inside x264 (ceiling about 2 %) |
| GNSS ±1 matched-filter banks on FPGA with shared chunks | 0 % (Yosys already merges them); the same bank as one matrix product on the M4: 3.2 times slower than FFT correlation |
| serial-commutator FFT instead of radix-2² single-path delay feedback | +0.8 % to +9.5 % slices, +56 % to 94 % flip-flops |
| streaming 8 × 8 SATD as delay-feedback stages; fully parallel 8 × 8 SATD | +15 % (ECP5), +24 % (iCE40); area × period -5.7 % |
| precision-planned radix-2^4 4,096-point FFT on ECP5-85F | area × period -9.2 %; the DSP saving is the known radix-2^4 design |
| two FFT streams per DSP48E2 at equal signal-to-noise | impossible: needs 2n + t + 1 ≥ 30 > 27 bits |
| filter-bank packing on 7-series DSP48E1 with whole-branch quantisers | alias edge -57.6 to -58.3 dB, needs -60 (reopened by graded groups, 6a) |
| one whole-branch packed filter-bank cascade reaching -80 dB far stop-band | the rounding floor alone is about -80.5 / -78 dB |
| global block floating point on a 16-bit FFT next to a strong tone | recovers 6-8 dB of the 30+ dB lost |
| shared-prefix maximum-likelihood decoders for short 3GPP codes | (32,8) and (20,K) below, (32,7) infeasible; the others about 0.15 % of a base station |
| shared-prefix decoder on iCE40 (no LUT RAM) | +86 % |
| the unregistered DSP-stage WHT on the UP5K at 48 MHz | about 46 MHz |
| ReLU-form absolute values in FPGA SATD | -2.8 % to +2.9 % |
| HEVC interpolator tricks as a motion-estimation LUT gain on Xilinx 7 | -7.6 % against a fair rival (6e) |
| before these waves: the 19/20-operation 8-point step inside FFHT's AVX2 code | 1 % to 24 % slower at every size up to 2^19 (longer dependency chain); published in wht-fewer-operations |
| before these waves: the low-operation-count circuits pipelined on iCE40 | the textbook design is smaller (363 against 677 cells at one adder per clock); published in wht-fewer-operations |

### 9. What was claimed and what survived

The later figure is the one used everywhere in this section.

| lead | first claimed | after the critics |
|---|---|---|
| matrix-unit FIR, fp32 | -75 / -87 / -68 % against vDSP and PFFFT | -59 / -53 / -50 % against Accelerate `sgemm` block-Toeplitz |
| matrix-unit 256-point FFT | -28 % | -28 % per core; whole chip -5 % to -12 % |
| fused (32,11) decoder | -35.0 % LUT4 | -31.3 % (stronger unfused baseline) |
| NEON decoder | -61 % | the new part is -18 % |
| Whisper spectrogram | -50.6 %; 4 cores -23.4 % (thread outputs unchecked) | -49.7 % to -50.9 %; 4 cores -22.9 % (checked); from files -1.6 % to -3 % |
| UP5K WHT | 246 against 633 cells at 60.5 MHz; 8k/16k points "a capability" | 253 against 633 at about 58 MHz; the larger sizes come from memory banks the baseline can use too |
| streaming 8 × 8 SATD | -35.9 % LUT4 | -24.2 % LUT4 against LUT-RAM transposition, +1.4 % slices |
| int16 rotate-and-quantise | -22.3 / -26.6 % | refuted (fp16 is as fast) |
| matrix-unit FIR bank | -36 / -41 % | -29 / -31 %, warm caches only |
| batched complex128 FFT, one core, in cache | -55.8 % (512), -61.9 % (4,096) | -41 % (512, against KFR), -44 % (4,096, against FFTW double with NEON, after the baseline rebuild; -52 % before it) |
| batched complex128 FFT, whole cluster | -28 / -31 % | -26 / -18 % against 4 × KFR |
| SEAL NTT, "default moduli only -8 %" | below the bar | BFV defaults at N = 2^12 to 2^14 -16.3 % to -20.8 % |
| telescope filter-bank packing | -50 % DSPs, as if new | the halving is WP486's; the new part saves about 17k LUTs per chip and 1.6-2.1 dB |
| streaming 4,096-point FFT on UP5K | 31.0 MHz, "no baseline fits" | 30.9 MHz over 12 seeds; rivals are ablations; inputs must stay inside the full-scale circle |
| optics step, one core | -56 / -49 / -62 % against PFFFT or KFR rows plus transpose | -38 / -47 / -62 % against FFTW's own 2-D plan |
| optics step, whole cluster | -19 / -16 / -21 % | tie at 1,024²; -12.7 % at 2,048²; -20.7 % at 4,096² |
| Welch, one core | -33 / -30 / -15 / -41 % | -33 / -27 / -11 / -41 % |
| Welch, four cores | -10.3 % (1,024), -11.2 % (4,096) | no size clears |
| Welch from Python | -83 % to -87 % against scipy and MNE | -62 / -66 % against a 6-line batched `rfft` Welch |
| carry-folded FFT | -23.9 % slices on NX; -52 % against a published core | about -20 % against the gated generic RTL; -52 % dropped; open tools only |
| folded UP5K FFT | -58 % (against a logic serial design), -82 % at equal rate; 61.3 MHz | -40 % against a 2-DSP design, -56 % against plain DSP inference; 50 MHz |
| HEVC interpolators as DSP chains | module -41.5 % LUT; motion estimation -20.4 % | module -18.6 % LUT, -18.2 % DSP; motion estimation -7.6 % LUT (below); later found to rest on an incorrect design (6e) |
| "no FFT on Arm's matrix extension exists" | prior-art claim | Arm's RAL library has merge requests adding 2-D FFTs with these outer products at 256 × 256 (their precision is not stated) |

### 10. Traps for anyone measuring on these chips and tools

- **SME code on the M4 needs `-mcpu=apple-m4`.** `-march=armv9*` builds crash with an illegal instruction.
- **FFTW double precision is scalar on Apple silicon unless configured with `--build/--host aarch64-apple-darwin`**
  (and `--enable-neon`), because `configure` checks for `aarch64` and Apple reports `arm64`. Without
  `--enable-armv8-cntvct-el0` FFTW has **no timer** on macOS/arm64, so its MEASURE and PATIENT plans are untuned
  (1.9 times slower at 1,024 points in one test). FFTW wisdom made on an efficiency core is bad wisdom (up to 72 %
  slower on a P-core); plan on a P-core.
- **PFFFT built by hand needs `-DPFFFT_ENABLE_NEON`**, or it is scalar.
- **Never use a shipped FFTW as an Arm baseline** (section 5).
- **On the M4, vDSP's FFTs share the matrix unit** with your own matrix-unit code: a "four-thread vDSP" baseline scales
  only 1.2 to 1.6 times, and a matrix-unit kernel shares the same limit.
- **Yosys's `synth_ice40 -dsp` silently rewrites hand-made SB_MAC16 DSP instances into plain multipliers.** But a
  logic-only rival with a plain `*` should go through `-dsp`, or it is unfairly large.
- **nextpnr times UP5K DSP paths with 0.1 ns placeholders.** Use [up5k_sta.py](programs/fpga-measure/up5k_sta.py),
  which re-times DSP and large-RAM paths from Lattice's characterised delays. **An SB_MAC16 multiply without its
  pipeline register is capped at 50 MHz by the datasheet**; turn on `PIPELINE_16x16_MULT_REG2`.
- **Gowin "slices touched" omits the slices that LUT-RAM blocks occupy.**
- **Gate-level simulation of CrossLink-NX netlists made with `synth_nexus -nobram`** gave wrong outputs for designs
  that are right in RTL, baseline and new design alike. Check the flow on the baseline first.
- **Yosys maps ECP5 MULT18X18D blocks with all their registers bypassed**, so pipeline registers stay in the fabric
  and inflate flip-flop counts. Yosys's CrossLink-NX flow does put them inside the DSP blocks; Lattice's own tools were
  not tried.
- **Background load matters.** An orphaned solver held one P-core through most of waves 11-14, so every "whole
  cluster" figure of those waves ran with three free P-cores, for both sides. Check what else runs before a cluster
  timing.

### 11. Wave 15 and the crash

Wave 15 started at about 05:40 on 12 October with eight developments: a whole-chip split-step optics loop, numpy's
convolve on the matrix unit (3d), the NTT inside OpenFHE on two chips, astropy's FFT convolution, the radar chain with
two real channels per FFT (6b), a 16,384-point vibration analyser on the UP5K (6c), a packed 4-bit correlator on ECP5,
and the HEVC motion estimation routed on ECP5 (6e). At about 06:35 the Mac ran out of memory and restarted. The cause
was ours: the convolve program's multi-thread path leaked its per-thread scratch memory on every call, a test made
thousands of calls, the disk was too full for swap to grow, and about 20 FPGA place-and-route runs were going at
the same time outside the job slots. The leak was fixed afterwards (each thread's scratch
now lives in a heap block freed when the thread ends; outputs bit-identical in all 3,471 cases tried), and a memory
guard was installed on the machine. Three reports were written (3d, 6c's 8,192-point FFT, 6e; the last still has two gaps) and one was half
written (6b). The OpenFHE, optics, astropy and correlator tracks left only working files, and no critic read any of
wave 15. No figure from the unfinished tracks is given here.

## What is still open, with a first experiment for each

1. **Why the SEAL port gains so much less on a non-Apple Arm core.** First experiment: time the new kernel and SEAL's
   own lazy forward NTT level by level on a Neoverse N2 (GitHub's Arm runners), and try a build tuned with
   `-mcpu=neoverse-n2`. If the lazy transform cannot gain there, use the new code on such cores only for the fully
   reduced transforms, where it clearly wins.
2. **The same NTT in OpenFHE, the more widely used library.** First experiment: route OpenFHE's NTT for primes below
   2^50 through the floating-point method, run OpenFHE's unit tests, and time its own CKKS benchmarks on the M4 and on
   an Arm server runner.
3. **A NEON build of FFTW for the people who install it from conda-forge or Homebrew.** First experiment: rebuild
   conda-forge's recipe with `--enable-neon` in both precisions and the host triple that turns double-precision NEON
   on, run FFTW's own `tests/check.pl` in both precisions, and time the probe's grid plus pyFFTW with the library
   swapped in. A 10-line recipe would help every Arm user of FFTW.
4. **The packed filter bank through AMD's tools.** First experiment: place and route one lane and one slice in
   Vivado at 200-250 MHz, and compare with AMD's own FIR Compiler as the rival.
5. **The carry-folded FFT through Lattice's tools.** First experiment: synthesise the folded and the gated generic
   designs with Lattice Radiant on CrossLink-NX. If Radiant already folds the multiplexers, the gain is gone.
6. **The convolve kernel, criticised and made safe.** First experiment: have a critic rebuild the baselines; put a C
   FFT overlap-save into the dispatcher above about 450 taps; add the per-cluster try-lock (3e); time 1, 4 and 8
   Python processes against numpy and scipy.
7. **The M4's second matrix unit is unused.** The efficiency cluster has its own, slower matrix unit, and no lead
   used it. First experiment: time the batched complex128 FFT on the efficiency cluster's unit alone, then one
   P-cluster thread plus one E-cluster thread plus NEON threads, against all-NEON on all ten cores.
8. **The optics step as a real loop.** First experiment: natural-order output inside the kernels and a 100-step
   split-step propagation loop, timed against `scipy.fft` with all workers and FFTW's threaded 2-D plan.
9. **pocketfft's batched vector path at large rows.** First experiment: cap the batched path's buffer at the L1 size
   and time 64 × 4,096 complex128 again; if it holds, it is a small patch to the FFT inside numpy and scipy.
10. **Power on small FPGAs.** No power figure exists anywhere in this section, though power decides sensor-class
    designs. First experiment: a switching-activity estimate with the open flow for the folded UP5K FFT and its
    rivals, then one measurement on a board.
11. **HEVC motion estimation routed at a usable clock, with the encoder's own SATD.** First experiment: pipeline the
    candidate-index path and the vertical filters (a latency-neutral split already exists and is cycle-exact), then
    route with 5 seeds.
12. **Finish wave 15's radar chain.** First experiment: place and route the packed and unpacked chains with 5 seeds
    each and write the area table.

## Programs and data

Every folder has a README with what it shows, how to build and run it, and what it needs. The programs keep their
original file headers, which name the round that wrote them; the Sources section maps those names to folders.
Third-party libraries are not included: FFTW (GPL) and KFR (GPL-2 or commercial) were used only for private
measurement, and PFFFT, pocketfft and FFHT are linked, not copied. Their build recipes are in the folder READMEs.

| folder | what it holds | re-run for this revision |
|---|---|---|
| [programs/README.md](programs/README.md) | what you need, and the build recipes of the third-party baselines (FFTW with NEON, PFFFT, pocketfft, FFHT, KFR) | - |
| [programs/timing-harness](programs/timing-harness/) | the shared timing harness `bench.h` and an example | used by every re-run below |
| [programs/baseline-wrappers](programs/baseline-wrappers/) | small C wrappers around KFR for the double-precision baselines | - |
| [programs/fpga-measure](programs/fpga-measure/) | `fpga_measure.py` (synthesis and place-and-route on 8 chips, harness, seeds), `fpga_equiv.py` (a formal check with SymbiYosys or a simulation of equality, with a broken copy), `up5k_sta.py` (UP5K DSP timing) | not re-run |
| [programs/neon-wht](programs/neon-wht/) | the NEON WHT baseline, its benchmark against FFHT, FFT baselines, vDSP examples | correctness check: identical |
| [programs/sme-fp32-fir](programs/sme-fp32-fir/) | the matrix-unit FIR (section 4) and the critic's `sgemm` baseline | correctness check: identical |
| [programs/sme-fp32-fft](programs/sme-fp32-fft/) | the matrix-unit 256- and 400-point FFTs (section 4) | correctness check: identical |
| [programs/sme-fp64-batched-fft](programs/sme-fp64-batched-fft/) | the batched complex128 FFT kernels (3a), the cluster and Python drivers, the baseline recheck | correctness check: identical; the broken copy fails |
| [programs/sme-optics-step](programs/sme-optics-step/) | the optics step (3b), the critic's FFTW 2-D comparison, the GPU check | not re-run (needs KFR) |
| [programs/sme-welch](programs/sme-welch/) | the one-pass Welch (3c), the critic's harness, the 6-line Python rival | not re-run |
| [programs/sme-cluster-trylock](programs/sme-cluster-trylock/) | the per-cluster try-lock (3e) with the log-mel and FFT workloads | not re-run |
| [programs/sme-fp64-convolve](programs/sme-fp64-convolve/) | the convolve drop-in (3d), with the memory-leak fix | correctness: 534 of 534 pass; broken copy fails 4 of 4; NaN pattern identical |
| [programs/fftw-neon-probe](programs/fftw-neon-probe/) | the shipped-FFTW probe and the pocketfft probes (section 5) | not re-run (needs the downloaded packages) |
| [programs/packed-pfb](programs/packed-pfb/) | coefficient searches, RTL generators, test benches, exhaustive unpacking checks (6a) | exhaustive 2^28 check: identical |
| [programs/carry-folded-sdf-fft](programs/carry-folded-sdf-fft/) | the FFT generator with folded carry cells, bit-true model, the critic's gated rival (6b) | carry-cell checks: 0 failures |
| [programs/up5k-folded-fft](programs/up5k-folded-fft/) | the folded 4,096-point and 25-bit 8,192-point UP5K FFTs, rivals, timing (6c) | precision model: identical |
| [programs/up5k-streaming-fft-and-wht](programs/up5k-streaming-fft-and-wht/) | the streaming 4,096-point FFT and WHT generators and models (6c) | FFT precision model: identical |
| [data/seal-ntt-port](data/seal-ntt-port/) | the SEAL benchmark tables on the M4 and on Neoverse N2, and the scripts that made them (section 2) | the tables are the original outputs |
| [data/packed-pfb](data/packed-pfb/) | the filter prototypes and the final quantised coefficient sets of 6a | the original files |

No timing was re-run for this revision: timings need a quiet machine and the timing lock, and the original numbers
with their critics' re-measurements are the evidence. The machine-readable summary is
[results.json](results.json).

## Credits

- **Microsoft SEAL** (Microsoft Research; MIT licence): the library the NTT result is measured in. **Intel HEXL**:
  the existing x86 hook that the port mirrors.
- **J. van der Hoeven and G. Lecerf**, "Implementing number theoretic transforms" (preprint, December 2024), and
  **FLINT**'s `fft_small`: the floating-point NTT arithmetic used in section 2.
- **FFTW** (Matteo Frigo and Steven G. Johnson; GPL), **KFR** (its authors; GPL-2 or commercial), **PFFFT** (Julien
  Pommier, with Hayati Ayguen and Dario Mambro in the marton78/pffft fork used here; BSD-like), **pocketfft** (Martin Reinecke; BSD-3),
  **FFHT** (Alexandr Andoni, Piotr Indyk, Thijs Laarhoven, Ilya Razenshteyn and Ludwig Schmidt; MIT; the copy inside
  PyTorch's ExecuTorch), **Apple Accelerate / vDSP** and **MLX**: the baselines.
- **conda-forge's `fftw-feedstock` issue #94, FFTW issues #251 and #129**, and the NixOS maintainers' measurements:
  the history behind section 5.
- **Yosys, nextpnr, icestorm, SymbiYosys** (YosysHQ and contributors), **Verilator** and **Icarus Verilog**: every
  FPGA number. The UP5K timing database is icestorm's, derived from Lattice's characterisation.
- **Xilinx WP486** (INT8 packing in DSP48E2): the known half of 6a. **G. Comoretto et al.**, the SKA-Low station
  firmware paper (arXiv 2002.10307): the channeliser's share of the chip. **Lattice patent US 9,390,210**: the prior
  art that may already cover 6b. **M. Garrido et al.**: the serial-commutator FFT tried in section 8. **nanamake's
  r22sdf** (MIT) and **ZipCPU's dblclockfft** (Dan Gisselquist): open FFT cores measured as references.
- **openasic-org's xk265** HEVC encoder: measured unchanged and never redistributed (6e). **x264**: the SATD test in
  section 8.
- **Arm's RAL library** (merge requests 44 and 46): prior SME2 outer-product FFTs. **KleidiAI** and **ExecuTorch**:
  the language-model kernels in section 7. **OpenAirInterface**: the base-station decoder named as the next rival.
- **GitHub's hosted Arm runners**: the Neoverse N2 measurements.

## Sources

Bookkeeping for readers who want to trace a number. The paths are those of my private working notes (not
published; the speed-wave folder is not in my research archive either).

- **Round notes** (12 October 2026, US Eastern): `wht42/LEADS-SPEED.md` (the running table of leads, last updated
  06:50 after wave 14; the best single source), `wht42/wave-11/NOTE.md` (00:56-01:03), `wave-12/NOTE.md`
  (02:42-02:50), `wave-13/NOTE.md` (04:30-05:00), `wave-14/NOTE.md` (06:20-06:40). Wave 15 has no note; its reports
  are `wave-15/av-1.md`, `av-4.md` (unfinished), `av-5.md`, `av-7.md`, and `av-picks.md` (the briefs).
- **Scouts** (01:03-01:48): `wht42/scout/scout-c0.md` (scientific Python shares), `scout-c1.md` (SEAL, big integers,
  language models), `scout-f0.md` (telescope F-engines, HEVC), `scout-f1.md` (SKA-Low channeliser, vendor FFT cores).
- **Baselines**: `wht42/baselines/README.md` (library builds), `baselines/FPGA.md` (the measuring tool),
  `baselines/STRONGEST.md` (the strongest processor FFT baselines, 04:40-05:25; the source of the -44 % correction in
  3a), `infra/inf-up5k-timing.md` (UP5K re-timing).
- **The SEAL port** (12-13 October; in the port's working folder, not published): `NOTES.md` (port, tests, review
  fixes), `bench/RESULTS.md` (M4 benchmark and control, about 02:20 on 13 October), the Arm-runner `RESULTS.md` (two
  workflow runs on 13 October 2026, 05:33-05:49 and 05:52-06:08 UTC, in a private repository), and the patch
  `neon-ntt-v4.6.2.patch` (10 files, +1,245 lines, sha256 `4eb1b5dc510ce3a3b1a00250a3b366bb7e805c59b0727b7fea66d81edd312d6c`;
  against SEAL tag v4.6.2, commit 31fc3194). The two skeptics' verdicts are recorded in my private notes of
  13 October, 02:30 US Eastern.
- **Agent labels to folders.** Developments, critics and fairness critics per wave: wave 11 aj-*, ak-cross, al-*;
  wave 12 am-*, an-cross, ao-*; wave 13 ap-*, aq-cross, ar-*; wave 14 as-*, at-cross, au-*; wave 15 av-*.
  - `programs/neon-wht`: aj-base-cpu. `programs/sme-fp32-fir`: aj-1, critic al-1. `programs/sme-fp32-fft`: aj-2
    (critic al-2). `programs/sme-fp64-batched-fft`: ap-0 (critic ar-0), baseline recheck inf-strongest.
    `programs/sme-optics-step`: as-1, critic au-0, fairness au-f. `programs/sme-welch`: as-0, critic au-1, fairness
    au-f. `programs/sme-cluster-trylock`: as-3, log-mel kernel ap-2. `programs/sme-fp64-convolve`: av-1.
    `programs/fftw-neon-probe`: au-h0 (probes P1-P3).
  - `programs/packed-pfb`: ap-4 (critic ar-3) and as-4. `programs/carry-folded-sdf-fft`: as-5, critic au-2.
    `programs/up5k-folded-fft`: as-6, critic au-3, and av-5. `programs/up5k-streaming-fft-and-wht`: ap-6 (critic
    ar-4) and am-5 (critic ao-1).
  - Section 2: ap-1 (critic ar-1, fairness ar-f), as-2, the port's round of 12-13 October (its exactness and build-shape skeptics, two rounds each). Section 6e: as-7,
    critic au-4, av-7. Section 7 and 8 entries: the corresponding rows of `LEADS-SPEED.md`.
- **Outside sources named above**: van der Hoeven and Lecerf, HAL 04841449, https://www.texmacs.org/joris/ntt/ntt.html;
  FLINT `src/machine_vectors.h`, https://github.com/flintlib/flint; conda-forge `fftw-feedstock` issue #94; FFTW
  issues #251 and #129 (https://github.com/FFTW/fftw3); Comoretto et al., arXiv 2002.10307; Xilinx WP486; US patent
  9,390,210; Arm RAL merge requests 44 and 46.
