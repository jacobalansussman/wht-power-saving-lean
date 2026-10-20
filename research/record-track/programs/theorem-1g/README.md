# Programs: Theorem 1G, the dilated lockstep tree, and its two skeptics

These programs check the construction of [../../proofs/theorem-1g.md](../../proofs/theorem-1g.md) in the job model.
They are Python 3.9 or newer, standard library only. Run each with `python3 -I -B` from this folder; each one took
under 40 seconds and under 300 MB on an Apple M4 Mac mini. All were re-run on 19 October 2026 for this revision, and
every output was identical to the one saved on 14 October (`../../data/theorem-1g/`).

The first four were written with the proof; the other five were written by two AI skeptics from the proof's text
alone, without reusing its code.

| program | what it does | command | output |
|---|---|---|---|
| `tree_flow.py` | builds the dilated lockstep tree for three test rules with grouped paid blocks, literally (explicit round sequences, L = 40 and 150) and merged (L up to 10^7); checks exact balance, sources in (L/m, L], Φ(exits) = Φ(sources), the exit bound and the unit bound | `python3 -I -B tree_flow.py` | `out_tree_flow.txt` |
| `recur.py` | the cost of the whole schedule (top split, tree, exits by fresh calls) on the test rules and the rank-1 control: T(n)/n^τ up to n | `python3 -I -B recur.py 800` and `python3 -I -B recur.py 6400` | `out_recur_800.txt`, `out_recur_6400.txt` |
| `recheck_t1.py` | re-runs the rank-1 τ-theorem's lockstep simulation and exact balance; imports `lockstep_sim.py` and `frontier_flow.py` from [../../../existence-routes/programs/tau-theorem-and-laser/](../../../existence-routes/programs/tau-theorem-and-laser/) | `python3 -I -B recheck_t1.py` | `out_recheck_t1.txt` |
| `hub_whatif.py` | the mixed price of the community's T74 and #1042 if every idle block were a head (the assumption later refuted by the head test), with the constants θ, η, κ of the theorem; control: the community's T74 figure | `python3 -I -B hub_whatif.py` | `out_hub_whatif.txt` |
| `indep_tree.py` | skeptic 1: an independent merged tree (toy rule and the T74 mixed rule, L up to 10^8) | `python3 -I -B indep_tree.py` | `out_indep_tree.txt` |
| `partial.py` | skeptic 1: the T74 price when only a fraction of the idle blocks can be cut and the rest are split 24 + 24 to stay separated; the smallest κ | `python3 -I -B partial.py` | `out_partial.txt` |
| `explicit_sim.py` | skeptic 2: explicit unit and job objects, random matching by size inside the prescribed classes, a generic topological sort of start, issue and finish events as a deadlock test; a "careless" matching across classes must deadlock | `python3 -I -B explicit_sim.py` | `out_explicit_sim.txt` |
| `cost_fit.py` | skeptic 2: its own cost recursion for the schedule; the third argument is the per-unit overhead as a multiple of W | `for r in D A B; do python3 -I -B cost_fit.py $r 12800 0,3; done` and `for r in C H J X; do python3 -I -B cost_fit.py $r 12800 0; done` | `out_cost_fit_12800_DAB.txt`, `out_cost_fit_12800_CHJX.txt` |
| `check_E.py` | skeptic 2: hypothesis (E) of Lemma A3 on the actual flows at every n (the last run, on the T74 rule up to 3,000,000, at every 997th n), the unit bound of Lemma E, the exit bound, and the ψ = 1 version (Lemma A4 predicts it fails at every n) | `for r in A B C H J; do python3 -I -B check_E.py $r 12800 1; done; python3 -I -B check_E.py T 20000 1; python3 -I -B check_E.py T 3000000 997` | `out_check_E.txt` |

What to expect: 0 unbalanced points and Φ ratio 1 everywhere; every lockstep run deadlock-free and every careless
run deadlocked; T(n)/n^τ flat for the rank-1 control and slowly growing, far below (ln n)^κ, for the grouped rules;
(E) never violated and its ψ = 1 version violated at every n. Rule X of `cost_fit.py` is not separated; it is
included as a first look at the open case (README of this section, section 6, item 2).

Status: these runs are *checked by computer* for the cases they cover; they do not prove the theorem, whose proof
is *argued by hand, not checked*.
