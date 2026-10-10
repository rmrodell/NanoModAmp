#!/usr/bin/env python3
"""Golden plan check V1: legacy counts on the full selection vs the paper counts tables.

Pass rule (docs/DECISIONS.md, G1-c 2026-10-06, refined after the first V1 run):
  exact      MPRA in vitro and in cellulo: totalReads, every *.count column and delrate equal.
  tolerance  endogenous: totalReads and every *.count column within max(3, 5% of paper totalReads);
             no separate delrate limit (one read at a ~26-read site moves delrate by ~0.04).
             Attributed causes: 20250418 = endoG1G2.fasta reference (F2; the original reference
             reproduces the paper exactly); 20251022 = unseeded umi_tools representative choice in the
             paper run (same dedup groups and depth, different kept read).
Every paper row for an included sample at a BED site must exist in the rerun, and vice versa.

Usage: v1_compare.py --v1 V1_DIR --out V1_results.tsv
"""
import argparse
import csv
import os
import sys

ENDO = "/scratch/users/rodell/NanoModAmp_build/endoBID_data"
POOL = "/oak/stanford/groups/nicolemm/rodell/BIDamplicon/Pool1"
# experiment: (paper table, rerun table relative to V1 dir, BED, rule)
EXPS = {
    "endo_20250418": (ENDO + "/20250418_BIDdetect_data.txt", "endo_20250418/counts/deldetect_factors.txt", "beds/endo_20250418.bed", "tolerance"),
    "endo_20251022": (ENDO + "/20251022_BIDdetect_data.txt", "endo_20251022/counts/BIDdetect_data.txt", "beds/endo_20251022.bed", "tolerance"),
    "mpra_invitro": (POOL + "/invitro/counts_delpos/BIDdetect_data_invitro_delpos.txt", "mpra_invitro/counts/BIDdetect_data.txt", "beds/mpra.bed", "exact"),
    "mpra_incell": (POOL + "/incell/counts_delpos/BIDdetect_data_incell_delpos.txt", "mpra_incell/counts/BIDdetect_data.txt", "beds/mpra.bed", "exact"),
}
META = ("celltype", "vector", "rep", "treat")
COUNTS = ("A.count", "C.count", "G.count", "T.count", "Deletion.count", "Insertion.count")


def load(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = {}
    for r in rows:
        key = tuple(r[m] for m in META) + (r["chr"], int(r["pos"]))
        n = int(float(r["totalReads"]))
        v = {"totalReads": n}
        for c in COUNTS:
            v[c] = int(float(r[c]))
        v["delrate"] = float(r["delrate"]) if r.get("delrate") not in (None, "", "NA") else v["Deletion.count"] / n
        out[key] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v1", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    out = open(a.out, "w")
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(("experiment", "rule", "celltype", "vector", "rep", "treat", "chr", "pos", "status",
                "paper_totalReads", "rerun_totalReads", "max_abs_count_diff", "abs_delrate_diff"))
    summary, all_ok = [], True
    for exp, (paper_f, rerun_f, bed_f, rule) in EXPS.items():
        sites = set()
        for line in open(os.path.join(a.v1, bed_f)):
            f = line.split("\t")
            sites.add((f[0], int(f[2])))  # BED end = 1-based position for both conventions in use
        rerun = load(os.path.join(a.v1, rerun_f))
        samples = {k[:4] for k in rerun}
        paper = {k: v for k, v in load(paper_f).items() if k[:4] in samples and (k[4], k[5]) in sites}
        n_ok = n_fail = 0
        worst_cnt = worst_dr = 0.0
        for k in sorted(set(paper) | set(rerun), key=str):
            p, r = paper.get(k), rerun.get(k)
            if p is None or r is None:
                status = "missing_in_rerun" if r is None else "missing_in_paper"
                w.writerow((exp, rule) + k + (status, p and p["totalReads"], r and r["totalReads"], "", ""))
                n_fail += 1
                continue
            dmax = max(abs(p[c] - r[c]) for c in COUNTS)
            ddr = abs(p["delrate"] - r["delrate"])
            worst_cnt, worst_dr = max(worst_cnt, dmax), max(worst_dr, ddr)
            if rule == "exact":
                ok = p["totalReads"] == r["totalReads"] and dmax == 0 and ddr < 1e-12
            else:
                lim = max(3, 0.05 * p["totalReads"])
                ok = abs(p["totalReads"] - r["totalReads"]) <= lim and dmax <= lim
            n_ok += ok
            n_fail += not ok
            w.writerow((exp, rule) + k + ("pass" if ok else "FAIL", p["totalReads"], r["totalReads"], dmax, "%.4g" % ddr))
        all_ok &= n_fail == 0
        summary.append("%-14s %-9s rows %3d  pass %3d  fail %3d  samples %2d  max|dcount| %g  max|ddelrate| %.4f"
                       % (exp, rule, n_ok + n_fail, n_ok, n_fail, len(samples), worst_cnt, worst_dr))
    out.close()
    print("\n".join(summary))
    print("V1 %s" % ("PASS" if all_ok else "FAIL"))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
