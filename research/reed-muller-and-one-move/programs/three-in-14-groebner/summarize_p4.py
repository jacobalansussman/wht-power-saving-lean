"""summarize_p4.py: coverage table of the four-array probe (probe4.py) from the fetched logs in jobs/p4/remote-out/.
A run stopped by hand has finished at least N candidates, N = its last "tried N" line (probe4.py prints it before
candidate N); for b5:0, whose copy printed no progress lines, N comes from p4_where.py (jobs/p4/where_b50_stop.txt)."""
import re, glob
R = "jobs/p4/remote-out/"
rows = [("b4 untargeted (|U| <= 6)", "ssh/log_b4.txt"), ("b4 targeted (|U| <= 7)", "ssh3/log_t_b4.txt"),
        ("b4plain untargeted (must find none)", "ssh3/log_b4plain.txt"), ("w3:0 untargeted (must find)", "ssh3/log_w3_0.txt")]
for k, f in [(0, "ssh3/log_t_b5_0.txt"), (1, "ssh/log_t_b5_1.txt"), (2, "ssh/log_t_b5_2.txt"), (3, "ssh/log_t_b5_3.txt"),
             (4, "ssh/log_c_b5_4.txt"), (5, "ssh2/log_t_b5_5.txt"), (6, "ssh3/log_c_b5_6.txt"), (7, "ssh2/log_t_b5_7.txt"),
             (8, "ssh2/log_c_b5_8.txt")]:
    rows.append((f"b5:{k} targeted (|U| <= 7)", f))
print("| run | candidates in the list | tried | filter passes (two orders UNIT) | result |")
print("|---|---|---|---|---|")
tot_list = tot_tried = 0
for name, f in rows:
    s = open(R + f).read()
    m = re.search(r"capped to the first (\d+) of (\d+)", s); n = int(re.search(r"candidate move sets (\d+)", s).group(1))
    full = int(m.group(2)) if m else n
    res = re.search(r"RESULT .*: (.*)", s)
    if res and res.group(1) != "None": tried = None; how = "until the first find (probe4.py stops there)"
    elif res: tried = n; how = f"{tried} (all{' of the cap' if m else ''})"
    elif "b5:0" in name:
        tried = int(re.search(r"matches \[(\d+)\]", open("jobs/p4/where_b50_stop.txt").read()).group(1)); how = f"{tried} (stopped)"
    else:
        tried = int(re.findall(r"tried (\d+) of", s)[-1]); how = f"at least {tried} (stopped)"
    passes = len(re.findall(r"^U=", s, re.M))
    print(f"| {name} | {full} | {how} | {passes} | {res.group(1) if res else 'stopped, nothing found'} |")
    if name.startswith("b5:"): tot_list += full; tot_tried += tried
print(f"\nblock 5, all nine sets: {tot_tried} of {tot_list} targeted candidates tried (at least), no filter pass")
