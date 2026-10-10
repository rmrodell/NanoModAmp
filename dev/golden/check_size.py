#!/usr/bin/env python3
"""Golden plan §4.6 step 2: measure the package size before V2.

Sums the capped FASTQs (select/<run>/capped), the references as they will ship (gzip -9),
the selected-chr BEDs, and a budget for expected outputs. Fails (exit 1) if FASTQs + references
exceed --data-cap-mb (default 22 MB, leaving 3 MB of the 25 MB cap for expected outputs).

Usage: check_size.py --golden GOLDEN_DIR [--report size_report.txt]
"""
import argparse
import glob
import gzip
import os
import sys

EXPS = {  # package experiment: (select runs, reference, BED)
    "endogenous": (["endo_20250418", "endo_20251022"], "/home/groups/nicolemm/rodell/fastas/endoG1G2.fasta", ["endo_20250418.bed", "endo_20251022.bed"]),
    "mpra_invitro": (["mpra_invitro"], "/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta", ["mpra.bed"]),
    "mpra_incell": (["mpra_incell"], "/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta", ["mpra.bed"]),
}
MB = 1e6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--golden", required=True)
    ap.add_argument("--report")
    ap.add_argument("--cap-mb", type=float, default=25.0)
    ap.add_argument("--data-cap-mb", type=float, default=22.0)
    ap.add_argument("--expected-budget-mb", type=float, default=3.0)
    ap.add_argument("--expect-fastqs", type=int, default=68, help="approved samples (G1-b); 0 = skip check")
    a = ap.parse_args()

    lines, data_total, bed_total, n_fq = [], 0.0, 0.0, 0
    lines.append("%-14s %6s %10s %10s %8s" % ("experiment", "fastqs", "fastq_MB", "ref_gz_MB", "bed_KB"))
    for exp, (runs, ref, beds) in EXPS.items():
        fqs = [f for r in runs for f in glob.glob(os.path.join(a.golden, "select", r, "capped", "*.fastq.gz"))]
        fq = sum(os.path.getsize(f) for f in fqs)
        ref_gz = len(gzip.compress(open(ref, "rb").read(), 9)) + os.path.getsize(ref + ".fai")
        bed = sum(os.path.getsize(os.path.join(a.golden, "v1", "beds", b)) for b in beds
                  if os.path.exists(os.path.join(a.golden, "v1", "beds", b)))
        data_total += fq + ref_gz
        n_fq += len(fqs)
        bed_total += bed
        lines.append("%-14s %6d %10.2f %10.3f %8.1f" % (exp, len(fqs), fq / MB, ref_gz / MB, bed / 1e3))
    total = data_total + bed_total + a.expected_budget_mb * MB
    ok = data_total <= a.data_cap_mb * MB and total <= a.cap_mb * MB
    if a.expect_fastqs and n_fq != a.expect_fastqs:
        lines.append("FASTQ count %d != expected %d" % (n_fq, a.expect_fastqs))
        ok = False
    lines.append("FASTQs + references: %.2f MB (limit %.1f MB)" % (data_total / MB, a.data_cap_mb))
    lines.append("+ BEDs + expected-output budget %.1f MB = %.2f MB (cap %.1f MB)" % (a.expected_budget_mb, total / MB, a.cap_mb))
    lines.append("SIZE %s" % ("PASS" if ok else "FAIL: lower the cap (nested) per G1-d"))
    text = "\n".join(lines)
    print(text)
    if a.report:
        open(a.report, "w").write(text + "\n")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
