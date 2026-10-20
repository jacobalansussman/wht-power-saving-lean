"""Write research/index.json, research/open-problems.md and the table of sections in research/README.md from the
sections' results.json files.

usage (from the root of the repository):  python3 -I -B research/make_index.py [--check]

Python 3.9 or newer, standard library only, under a second. With --check it writes nothing and exits 1 if either
file differs from what it would write (so the two files can be checked against the sections).
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORDER = [
    "moonshot-small-blocks", "moonshot-block4-block5", "reed-muller-and-one-move", "frame-packing-two-arrays",
    "why-the-savers-save", "block-bounds-and-new-rules", "theory-of-theories-and-logic", "existence-routes",
    "other-fields", "idea-rounds", "speed-waves", "record-track",
]
WORDS = ["proved in Lean", "checked by computer", "argued by hand, not checked", "refuted", "open"]
SUMMARY = {
    "moonshot-small-blocks": "The moonshot waves at blocks 2 and 3: no unit in the everything-free game or among mixed "
                             "units; block 2 never saves in any game with bounded free steps (Theorem 2W); every unit "
                             "with two or more arrays needs at least W + m moves.",
    "moonshot-block4-block5": "The moonshot waves at blocks 4 and 5: no unit; the 2-in-7 near miss at block 4 is still "
                              "open by 640 out of 137,280; many shapes excluded under stated hypotheses.",
    "reed-muller-and-one-move": "The Reed-Muller family (whose saving would tend to 1) and units that save a single "
                                "move: neither found; where the family is dead, and what is still undecided.",
    "frame-packing-two-arrays": "How a unit saves at all: every saving is frame packing; the frame game holds the whole "
                                "problem; whether any two-array unit saves (Proposition 6) is open.",
    "why-the-savers-save": "Where the record circuits' saving comes from (one free direction per array, a shared centre), "
                           "and why every lever tried on that mechanism failed to lift the record.",
    "block-bounds-and-new-rules": "Lower bounds on the moves a unit needs, three of them proved in Lean (Theorem R, block "
                                  "2, block 3), the ceiling of the method, and which new free steps (Toffoli, column "
                                  "mixing) are allowed and what they buy.",
    "theory-of-theories-and-logic": "Which kinds of lower-bound proof and construction can exist at all: the growth "
                                    "exponent screen, the spectral and angle factors, and what logic says about "
                                    "existence proofs.",
    "existence-routes": "Proving that a unit exists without finding it: a unit known only to exist counts (proved in "
                        "Lean), but no route reached an open shape; which routes work as certifiers and which are blind.",
    "other-fields": "How the problem relates to famous problems, quantum circuits and network coding: exact "
                    "reshapings, the textbook dictionary, and what each field's tools prove here.",
    "idea-rounds": "The two thinking rounds (10 October evening, 16 October): closed doors, structural facts and the "
                   "leads still open.",
    "speed-waves": "Real processors and FPGAs (11-13 October): one application-level gain (an NTT inside Microsoft SEAL, "
                   "21.2 % faster), kernel gains, and every lead that failed.",
    "record-track": "The race for the record (13-20 October): the Lean ladder to 1.6257e-3, the levers that died, the "
                    "recursion lever (Theorem 1G, heads, the half cut) and the frontier.",
}


def fmt_worth(w):
    if w is None:
        return ""
    if isinstance(w, (int, float)):
        return f"{w:.3g}" if w >= 0.1 else f"{w:.2e}".replace("e-0", "e-")
    if isinstance(w, dict):
        return "; ".join(f"{k}: {fmt_worth(v)}" for k, v in w.items())
    return str(w)


def cell(s):
    return " ".join(str(s).split()).replace("|", "\\|")


def build():
    sections, problems = [], []
    for key in ORDER:
        d = json.loads((HERE / key / "results.json").read_text())
        counts = {w: 0 for w in WORDS}
        pending = 0
        for f in d.get("findings", []):
            s = f.get("status")
            if s in counts:
                counts[s] += 1
            else:
                pending += 1
        entry = {
            "key": key,
            "title": d.get("title", key),
            "path": f"research/{key}/",
            "dates": d.get("dates", ""),
            "summary": SUMMARY[key],
            "findings": len(d.get("findings", [])),
            "counts": counts,
            "open_questions": len(d.get("open", [])),
            "programs": len(d.get("programs", [])),
            "results": f"research/{key}/results.json",
            "readme": f"research/{key}/README.md",
        }
        if pending:
            entry["findings_with_status_pending"] = pending
        sections.append(entry)
        for o in d.get("open", []):
            problems.append((key, o))
    index = {
        "schema": "wht-research-index/1",
        "repository": "jacobalansussman/wht-power-saving-lean",
        "revision": "fifth revision",
        "covers": "2026-10-10 to 2026-10-20",
        "made_with": "AI agents (Claude) directed by the author",
        "human_review": "none: no person has read any of the proofs",
        "status_words": WORDS,
        "related": [{"name": "wht-fewer-operations",
                     "url": "https://github.com/jacobalansussman/wht-fewer-operations",
                     "holds": "the 5/6 circuit result and the ten idea waves of 10-11 October 2026"}],
        "open_problems": "research/open-problems.md",
        "generated_by": "research/make_index.py",
        "sections": sections,
        "totals": {
            "findings": sum(s["findings"] for s in sections),
            "counts": {w: sum(s["counts"][w] for s in sections) for w in WORDS},
            "open_questions": len(problems),
        },
    }
    lines = [
        "# Open problems of the research since the fourth revision",
        "",
        "Every open question of the twelve sections, each with the first experiment its section proposes. This file is",
        "generated from the sections' `results.json` by `make_index.py`; the sections' READMEs give the context, the",
        "hypotheses and the status of what is known around each question. \"Worth\" is the saving the answer would give",
        "if it went the hopeful way (a price, not a proved figure), where the section states one.",
        "",
    ]
    for s in sections:
        qs = [o for k, o in problems if k == s["key"]]
        if not qs:
            continue
        lines += [f"## [{s['title']}]({s['key']}/README.md)", "",
                  "| | question | worth | first experiment |", "|---|---|---|---|"]
        for i, o in enumerate(qs, 1):
            lines.append(f"| {i} | {cell(o.get('question', ''))} | {cell(fmt_worth(o.get('worth', o.get('worth_if_it_existed'))))} "
                         f"| {cell(o.get('first_experiment', ''))} |")
        lines.append("")
    dates = {"moonshot-small-blocks": "14-15 Oct", "moonshot-block4-block5": "14-16 Oct", "reed-muller-and-one-move": "15-16 Oct",
             "frame-packing-two-arrays": "15 Oct", "why-the-savers-save": "13-15 Oct", "block-bounds-and-new-rules": "13-16 Oct",
             "theory-of-theories-and-logic": "14-15 Oct", "existence-routes": "14-15 Oct", "other-fields": "10-15 Oct",
             "idea-rounds": "10 and 16 Oct", "speed-waves": "11-13 Oct", "record-track": "13-20 Oct"}
    about = {"moonshot-small-blocks": "the moonshot waves at blocks 2 and 3: the everything-free game and mixed units",
             "moonshot-block4-block5": "the moonshot waves at blocks 4 and 5: the 2-in-7 near miss and the shape searches",
             "reed-muller-and-one-move": "the Reed-Muller family, and units that save a single move",
             "frame-packing-two-arrays": "how savings happen: crossing slack, frame packing, the frame game, Proposition 6",
             "why-the-savers-save": "why the record circuits save, and what caps them: anatomy, centre, copy wall, helpers",
             "block-bounds-and-new-rules": "lower bounds at each block size, in Lean; Toffoli and column-mixing maps",
             "theory-of-theories-and-logic": "what our methods can and cannot do: growth exponents, certificates, logic",
             "existence-routes": "proving that a unit exists without finding it",
             "other-fields": "famous problems, network coding and quantum computing",
             "idea-rounds": "the two thinking rounds: closed doors, structural facts, open leads",
             "speed-waves": "speed-first waves 11 to 15: real processors and FPGAs",
             "record-track": "the Lean ladder, the levers that failed, the recursion lever and the frontier"}
    rows = ["| section | dates (2026) | what it is about | findings: proved in Lean / checked by computer / argued by hand / "
            "refuted / open | open questions |", "|---|---|---|---|---|"]
    for s in sections:
        c = " / ".join(str(s["counts"][w]) for w in WORDS)
        if s.get("findings_with_status_pending"):
            c += f", and {s['findings_with_status_pending']} whose Lean status is given with the result"
        rows.append(f"| [{s['key']}]({s['key']}/README.md) | {dates[s['key']]} | {about[s['key']]} | {c} | {s['open_questions']} |")
    return json.dumps(index, indent=1, ensure_ascii=False) + "\n", "\n".join(lines), "\n".join(rows)


def readme_with(table):
    p = HERE / "README.md"
    t = p.read_text()
    a = t.index("<!-- the table below is written by make_index.py -->\n") + len("<!-- the table below is written by make_index.py -->\n")
    b = t.index("\n<!-- end of the table written by make_index.py -->")
    return t[:a] + table + t[b:]


if __name__ == "__main__":
    js, md, table = build()
    targets = [(HERE / "index.json", js), (HERE / "open-problems.md", md), (HERE / "README.md", readme_with(table))]
    if "--check" in sys.argv:
        bad = [p.name for p, t in targets if not p.exists() or p.read_text() != t]
        print("differs: " + ", ".join(bad) if bad else "index.json, open-problems.md and the README table are up to date")
        sys.exit(1 if bad else 0)
    for p, t in targets:
        p.write_text(t)
    print(f"wrote index.json ({len(json.loads(js)['sections'])} sections), the README table and open-problems.md "
          f"({json.loads(js)['totals']['open_questions']} questions)")
