#!/usr/bin/env python3
"""Build the single endogenous BIDdetect_data.txt that Figure 2 analysis_endo.R reads.

No combined Figure 2 table was found on Oak, so the harness constructs one:
  - one input table per run (formatted: celltype vector rep treat chr pos ... [delrate]);
  - vector labels normalised to WT/KD (Becca 2026-10-06: P36/pLKO = WT, P97/shPUS7 = KD;
    20251022 already uses WT/KD) — analysis_endo.R filters on vector == "WT"/"KD";
  - rep made unique per run as <run>_<rep>, since rep k is batch-paired within a run only (D17);
  - delrate = Deletion.count / totalReads where the table lacks it (20250418 DelDetect output);
  - optional filters to the selected chrs and samples (used for the paper-table test).
Output columns follow the 20251022 table (with delrate).

Usage: legacy_combine_endo.py OUT RUN=TABLE [RUN=TABLE ...] [--chrs FILE] [--keep-reps RUN:1,3 ...]
"""
import argparse
import csv

VEC = {"P36": "WT", "pLKO": "WT", "P97": "KD", "shPUS7": "KD", "WT": "WT", "KD": "KD"}
COLS = ["celltype", "vector", "rep", "treat", "chr", "pos", "gene", "totalReads", "A.count", "C.count",
        "G.count", "T.count", "Deletion.count", "Insertion.count", "ref", "kmer", "strand", "delrate"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("tables", nargs="+", help="RUN=path")
    ap.add_argument("--chrs")
    ap.add_argument("--keep-reps", nargs="*", default=[], help="RUN:1,3 keeps only those reps of RUN")
    a = ap.parse_args()
    chrs = {l.strip() for l in open(a.chrs) if l.strip()} if a.chrs else None
    keep = {k.split(":")[0]: set(k.split(":")[1].split(",")) for k in a.keep_reps}
    rows = []
    for spec in a.tables:
        run, path = spec.split("=", 1)
        for r in csv.DictReader(open(path), delimiter="\t"):
            if chrs is not None and r["chr"] not in chrs:
                continue
            if run in keep and r["rep"] not in keep[run]:
                continue
            r["vector"] = VEC[r["vector"]]
            r["rep"] = "%s_%s" % (run, r["rep"])
            if not r.get("delrate"):
                r["delrate"] = repr(int(r["Deletion.count"]) / int(r["totalReads"]))
            rows.append([r[c] for c in COLS])
    with open(a.out, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        w.writerows(rows)
    print("wrote %d rows to %s" % (len(rows), a.out))


if __name__ == "__main__":
    main()
