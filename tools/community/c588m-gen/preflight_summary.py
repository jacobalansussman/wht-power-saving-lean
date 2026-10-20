#!/usr/bin/env python3
"""m-gen: one table of every exact Python preflight of the C588M build, each with its verdict (read from the logs).
usage: preflight_summary.py <m-gen dir>"""
import sys, os, json

M = sys.argv[1]
L = lambda f: json.load(open(os.path.join(M, "logs", f)))
rows = []


def add(name, ok, detail):
    rows.append(dict(check=name, pass_=bool(ok), detail=detail))


p = L("prep_cert.json")
add("input certificate = Shea's candidate with the 112 copies folded back = x-mutable's file; histogram replayed",
    p["equal_to_x_mutable"] and p["N"] == 5916 and p["R"] == 504,
    dict(sha256=p["out_sha256"], R=p["R"], N=p["N"], hist=p["hist"], entrance_rank_sum=p["entrance_rank_sum"]))
a = L("scalar_audit_z.json")
add("rational scalar audit, mutable kinds (x: hid, y: identity, e: E7); Shea's rule gives SPar",
    a["SPar_by_Shea_rule"] == [10, 21, 4, 8, 8] and all(m["pass"] for m in a["checks"]),
    dict(SPar=a["SPar_by_Shea_rule"], largest_digit=a["largest_uncancelled_digit"],
         x_rows_restored=a["checks"][0]["x_rows_restored_(classic_hx)"]))
z = L("zmirror-c588m.json")
add("exact packed mirror of the kernel terms: scalCheckZ x4 (10,21,4,8,8), yCheckZ, eCheckZ x3 (10,8,4,8,8) all true",
    z["main"]["all_true"], dict(scal=[x["value"] for x in z["main"]["scal"]], ychk=z["main"]["ychk"]["value"],
                                echk=[x["value"] for x in z["main"]["echk"]]))
add("margin: at sw - 1 the scalar batches fail (digit width is tight, as Shea's 21)",
    not any(x["value"] for x in z["sw_minus_1"]["scal"]), dict(scal_sw20=[x["value"] for x in z["sw_minus_1"]["scal"]]))
c = L("mirror-controls.json")
add("mirror controls vs facts Lean decided in Shea's build (classic accept at sw 21, his bad scalar false, sw 20 false); "
    "copy-free cert fails the classic checks",
    c["shea_classic_10_21_1_8_8_batch0"]["value"] and not c["shea_bad_scalar_rejected_expect_false"]["value"]
    and not c["shea_classic_sw20_batch0"]["value"] and not c["ours_classic_expect_false"]["all_true"],
    dict(ours_classic_first_failure=c["ours_classic_expect_false"]["scal"][0].get("why")))
g = open(os.path.join(M, "logs/generate_local.log")).read()
add("Sussman's Python mirror inside gxgen: label side E0/E1/E2/E6 (frames, pairs, replay, scatter cut, histogram) and "
    "the mutable-kind scalar side accept", "C588M: mirror accepts" in g,
    [l for l in g.splitlines() if "mirror accepts" in l])
r = L("roundtrip.json")
add("generated Lean data parse back to the certificate; kernel-term mirror on the parsed text all true",
    r["parsed_equals_certificate"] and r["kernel_terms_on_parsed_text"]["all_true"], r["kernel_terms_on_parsed_text"])
for v in ("annealed", "cyclic", "packed"):
    t = L("trie_check-%s.json" % v)
    add("tries of the generated data (%s): tab codes, st0 = start, res_ok (schedRes), cen_ok" % v, all(t.values()), t)
    pr = json.load(open(os.path.join(M, "out/price/price-%s.json" % v)))
    add("bank table %s: schedK mirror true; duplicate owner and wrong pad controls false; deficit 2560" % v,
        pr["schedK_mirror"] and not pr["schedK_duplicate_owner_control"] and not pr["schedK_wrong_pad_control"]
        and pr["deficit"] == 2560,
        dict(banks=pr["banks"], W=pr["W"], pads=pr["pad_hist"]))
    gl = json.load(open(os.path.join(M, "out/global-%s/logs/global-C588M.json" % v)))
    add("rate certificate %s (Sussman's foldlib.plan_block, B = 10, s = 40) at the kernel numerator; DFT = kernel - 1" % v,
        gl["kernel_numerator"] == pr["kernel_numerator"] and gl["dft_numerator"] == pr["kernel_numerator"] - 1,
        dict(root_bracket=[pr["exact_rate"]["root_lower_decimal"][:18], pr["exact_rate"]["root_upper_decimal"][:18]],
             kernel=gl["kernel_numerator"], dft=gl["dft_numerator"], rational_slack=gl["rational_slack"]))
k = L("control_facts.json")
add("Lean negative controls predicted false (scatter sign flipped; classic check; duplicate owner), positives true",
    (not k["bad_scalar_rejected"]["value"]) and (not k["classic_rejected"]["value"]) and k["positive_scal0"]["value"]
    and k["positive_sched"] and not k["duplicate_owner_rejected"],
    dict(note="Shea's own control (first A gate doubled) would be TRUE here, so it is replaced",
         shea_style=k["shea_style_first_A_gate_doubled_NOT_rejected"]["value"]))
x = L("relocate_xmut.json")
add("x-mutable's 7 engine modules relocated to Work/XMut (import prefix only)", len(x) == 7,
    {f: v["imports_rewritten"] for f, v in x.items()})
out = dict(all_pass=all(r["pass_"] for r in rows), checks=rows)
print(json.dumps(out, indent=1))
