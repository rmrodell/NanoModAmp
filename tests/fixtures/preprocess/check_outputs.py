#!/usr/bin/env python3
"""Check harness outputs against expected.json (WP2 tests 1-11, 13, 14 on the fixtures).

Usage: check_outputs.py <outdir> <sample> [<sample> ...]
Exit 0 if every check passes; prints one line per check.
"""
import csv
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
exp = json.load(open(os.path.join(HERE, "expected.json")))
outdir, samples = sys.argv[1], sys.argv[2:]
fails = 0


def check(name, ok, detail=""):
    global fails
    print("%s  %s%s" % ("PASS" if ok else "FAIL", name, ("  (" + detail + ")") if detail and not ok else ""))
    fails += not ok


def reads(sample):
    with gzip.open(os.path.join(outdir, "reads", sample + ".fastq.gz"), "rt") as fh:
        lines = fh.read().splitlines()
    return [(lines[i][1:].split()[0], lines[i + 1]) for i in range(0, len(lines), 4)]


def funnel(sample):
    with open(os.path.join(outdir, "funnel", sample + ".read_funnel.tsv")) as fh:
        return {r["step"]: int(float(r["records"])) for r in csv.DictReader(fh, delimiter="\t")}


O1, E1 = exp["inserts"]["O1"], exp["inserts"]["E1_INS"]
u = exp["umis"]
for s in samples:
    if s == "S_dir":
        names = [n for n, _ in reads(s)]
        check("S_dir: directory concatenated in lexical order (test 8)",
              [n.split("_")[0] + "_" + n.split("_")[1] for n in names] == exp["dir_order"], str(names))
        continue
    f = funnel(s)
    want = exp["funnel"][s]
    got = {k: v for k, v in f.items() if k != "mapped"}
    check("%s: funnel steps and records (tests 10, 11, 13, 14)" % s, got == want, "got %s want %s" % (got, want))
    check("%s: mapped >= primary (test 7)" % s, f.get("mapped", 0) >= f.get("primary", 0))
    if s in exp["both_orientations"]:
        with open(os.path.join(outdir, "both", s + ".tsv")) as fh:
            row = next(csv.DictReader(fh, delimiter="\t"))
        check("%s: reads in both orientations" % s,
              int(row["n_reads_in_both_orientations"]) == exp["both_orientations"][s], str(row))
    r = dict(reads(s))
    if s == "mpra_ont":
        sense = [seq for n, seq in r.items() if n.startswith("sense_u2")]
        anti = [seq for n, seq in r.items() if n.startswith("antisense_u3")]
        check("mpra_ont: sense read trimmed to exactly the insert (test 1)", sense == [O1], str(sense))
        check("mpra_ont: antisense read identical to its sense twin after RC (test 2)", anti == [O1], str(anti))
        dropped = [n for n in r if n.split("_")[0] in ("only", "too")]
        check("mpra_ont: one-adapter and too-short reads discarded (tests 3, 4)", not dropped, str(dropped))
        check("mpra_ont: UMI appended to the read name with '_' (test 5)",
              any(n == "sense_u2_" + u["u2"] for n in r), str(sorted(r)))
    if s == "endo_ont":
        check("endo_ont: sense read trimmed to exactly the insert, UMI removed (test 1)",
              r.get("e_u2_" + u["uE2"]) == E1, str(r.get("e_u2_" + u["uE2"])))
    if s == "mpra_pool":
        check("mpra_pool: pool-flanked reads kept, no UMI in names (test 13)",
              sorted(n for n in r) == ["p_antisense", "p_sense1", "p_sense2"] and all(seq == O1 for seq in r.values()),
              str(sorted(r)))
    if s == "endo_3prime":
        names = sorted(n.rsplit("_", 1)[0] for n in r)
        check("endo_3prime: 3′-only and adapter-less reads kept, short read dropped (test 14)",
              names == ["t_3p_u1", "t_3p_u2", "t_no_adapter"], str(names))
        check("endo_3prime: 3′ adapter trimmed, UMI extracted (test 14)",
              r.get("t_3p_u1_" + u["uE1"]) == E1, str(r.get("t_3p_u1_" + u["uE1"])))
print("%d check(s) failed" % fails if fails else "all checks passed")
sys.exit(1 if fails else 0)
