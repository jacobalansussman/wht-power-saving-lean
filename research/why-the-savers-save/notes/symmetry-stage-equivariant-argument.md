> Working note written by an AI agent during the round of 13 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# d-I3-stage-equivariant: can a total cross a stage boundary without a copy?

13 October 2026. Standard library only. Every number below is exact.

## Answer

No, not on the hub's complex nine-stage geometry, and not for any decoder. The obstruction comes from the
routing, not from the decoder, so no symmetry of the decoder can remove it.

- The boundary maps A_T come exactly from the hub's own files.
  - The routes rho_{s,q} are rebuilt from the hub's written specification
    (`split-centres-global/proofs/INHERITED-COMPLEX-GEOMETRY.txt`) and from the chart code in the hub's checkers
    (`punctured112-global/code/audit_112_independent.py`; the same lines are in #1042's `t16/indep16.py`).
  - They match the hub's `DATA-PHASES.json` on all 120 ports: the idle-child axes on all 120, and all 1,080 stage
    translations (exit = e_z + rho_s(e_z)).
  - The #1042 k8 word uses the same `DATA-PHASES.json` and the same charts.
- **What the boundary does.** Invocation d at stage s processes port q's piece of data at class d·rho_{s,q}. The
  template's route semantics say: "data: physical cover class d*rho(stage,port)"; "work: private family at class d,
  erased before reuse". So the stage-(s+1) invocation that continues that piece is d·A_q, where
  A_q = rho_{s,q}·rho_{s+1,q}.
- **Transportable totals.**
  - Hub decoder: 0 of 8, at each of the 5 hand-over boundaries. Each total's 56 ports continue into 56 different
    invocations.
  - Best over the decoder's orbit: 0. This needs no search, because A_q depends only on the port's label.

## The table of A_T (full matrices in `out/A_table.json`)

| boundary | A_q | distinct over the 112 ports | 9×9 local block | A_q(E_{s+1}) | order | coordinates moved |
|---|---|---|---|---|---|---|
| 0→1 | rho_1 | 112 | q qᵀ | q^⊥ + E_0 | 2 | 17 |
| 1→2 | rho_1 rho_2 | 112 | q qᵀ | q^⊥ + E_1 | 3 | 25 |
| 3→4 | rho_3 rho_4 | 112 | q qᵀ | q^⊥ + E_3 | 3 | 25 |
| 5→6 | rho_5 rho_6 | 112 | q qᵀ | q^⊥ + E_5 | 3 | 25 |
| 7→8 | rho_7 rho_8 | 112 | q qᵀ | q^⊥ + E_7 | 3 | 25 |

Checks at every boundary and every port:
- The boundary identities hold: target end = next zero frame, and full frame = next source line. In other words,
  the data levels are X_s = U + E_s and Y_s = q^⊥ + E_s. All 8 × 112 cases pass.
- The scalar hand-over identity holds at 5 boundaries: 0→1, 1→2, 3→4, 5→6 and 7→8.
  - The content is exact at 0→1 and 7→8.
  - At 1→2, 3→4 and 5→6 a bridge term comes in between.
  - This is a replay of the template's own recipe through d-g-crossstage's `stagemap.py`.

## The obstruction (proved, and route-independent)

**Claim.** Take any choice of routes g_{s,T} that keeps a data piece's frames continuous across the boundary, that
is, with no paid move between the two stages. Then the continuation factor A_T = g_{s,T}·g_{s+1,T}⁻¹ satisfies

A_T(E_{s+1}) = q_T^⊥ + E_s.

**Proof.** A piece's frame in its own class coordinates is g⁻¹(window frame). The target-end identity says
g_{s,T}⁻¹(q_T^⊥ + E_s) = g_{s+1,T}⁻¹(E_{s+1}). Apply g_{s,T} to both sides.

**Consequences.**
1. The subspace q^⊥ + E_s determines q: intersect it with U and take the perpendicular. Checked on all 120
   labels.
2. Right multiplication acts freely on the cover. So two ports with different labels never continue into the same
   invocation, whatever routes are chosen.
3. A register built in one invocation from two or more ports with different labels never has all its pieces
   together in any single later invocation. This covers a total, and also a target pre-loaded with a total.
4. On P112 the 112 labels are distinct. So the only transportable "totals" are those that read a single port,
   which are just the data registers the hub already carries across. A decoder made only of such totals costs at
   least 112.

This proves the narrow lemma of `../../copywall/NOTE.md`, using only the hub's own route semantics and no new rule:
"a total built in stage i cannot be reused through stage i's addressing".

The brief's step 2 test was "restrictions agree modulo the Hadamard-free Clifford group". As stated it does not
decide transportability:
- linear address maps are free, so any two injective restrictions agree modulo that group;
- the deciding invariant is the cover class d·A_q, which is exact.

## The open escape from the copy-wall note, closed

The escape was: "pre-load the sum into a target register whose climb is already paid".
- **Frames.** It needs P_a ⊆ q_S^⊥, which never holds: P_a = q_a^⊥ for a punctured label q_a. That is 0 of 896
  (total, host) pairs.
- **Data.** In the host's continuation invocation, at most 1 of the pre-loaded pieces is present, and that one is the
  host's own. The other pieces are foreign data, and they would have to be cancelled later.
- **Frame test from the stage map.** For 0 of 448 (total, port) pairs does a total's frame reach the new window.
- **Paid re-routing.** Moving the pieces so that one total's 56 ports share a continuation costs at least 220 moves
  per total per boundary, against the 8 its star copy costs. This uses the generous cost for steps between frames
  that are not nested. Under the hub's strict rule that frames only climb, such moves are not allowed at all.

## The small case, settled exactly first (`small_case.py`)

- **(A)** h = 3, m = 5. Exhaustive over all 720² pairs of routes in O(5,2):
  - each port has 2,880 admissible pairs;
  - every admissible A maps E_1 onto q^⊥;
  - the sets of A for different ports overlap in 0 elements.
- **(B)** h = 4, m = 7, exhaustive over O(7,2) (1,451,520 elements) with g_0 = I. NOT RUN TO THE END: the pure-Python
  enumeration ran 13 minutes under a memory-capping wrapper and I stopped it at the time box. The command is still in
  `small_case.py B`. It adds nothing that (A) and the two-line proof do not already settle.
- **(C)** h = 3, m = 7, using the hub's routes, with the cover as the invocations. This is a symbolic replay of a total
  over two ports pre-loaded into a third port's target, at stage 1 → 2. The run covers d = I and 2,000 random
  invocations of O(7,2) (`out/small_case_Cs.txt`). In 0 of 6,003 (invocation, host) cases does the host's
  continuation hold a piece the total carries.
  - Left multiplication by d commutes with every route, so d = I alone already decides the question.
  - The whole-cover run is `small_case.py C`. It is slow and was not completed.

## Control

`pricer.py` reproduces #1042 exactly: b = 0.001674989015937672328468, kappa = 0.001672188119205421, certified
(`out/control1042.txt`). Nothing was priced as a gain, because no legal route exists.

## Reproduce

```sh
python3 -I stage_equivariant.py          # 3 s; writes out/stage_equivariant.txt, out/summary.json, out/A_table.json
python3 -I small_case.py A Cs            # 22 s
python3 -I small_case.py B C            # whole O(7,2): over 13 min, not completed (run under a memory cap)
python3 -I ../../bound-check/m-anatomy/pricer.py --profile '{"1":349577,"2":60193,"3":62608,"4":20448,"5":11808,"7":3,"8":1152,"16":47,"48":4850}' --m 73 --W 14336 --deficit 5120
```

## Scope

- The proof is about the hub's compressed-window cover construction, for any routes that keep frames continuous.
  It is not about the frame-game rules for additions, so it also holds under our Lean engine's wider rules,
  provided the same global cover is used.
- A different global construction, one not built from a cover of charts with one window per invocation, is outside
  its scope.

## Credit

- The stage map and hand-over identity come from our copywall d-g-crossstage.
- The routes, charts, route semantics and data phases are the hub's:
  - PR514's complex-nine-native geometry;
  - icekylinx's punctured112-global independent audit;
  - Hao Yan's (@FanPUYen) #858 T16 checker in #1042.
- The decoder and the hand-over identity also come from the hub's format-limits-corrected notes, #446,
  DreamingOfClouds #717 and the #799 lineage.
- #1042 is the T16 owner packing on #1018's seven-slot k8 word.
