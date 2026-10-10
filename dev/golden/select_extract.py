#!/usr/bin/env python3
"""Golden plan §4.3–4.5 for one sample: read-ID assignment, background sampling, extraction.

Read IDs come from the rerun pre-dedup BAM (the legacy scripts apply no primary/MAPQ
filter, F1): every read with any alignment on a selected chr is selected, with all its PCR
duplicates. Each read is assigned to the chr of its primary alignment (or its first
selected-chr record if the primary is elsewhere).

Two FASTQs are written, both raw records in original order with adapters intact:
  full/<sample>.fastq.gz     all selected reads + background   (input to V1)
  capped/<sample>.fastq.gz   ≤ cap selected reads per chr + background (the package, G1-d)
Sampling is seeded and nested: the capped set for cap N is a prefix of a fixed per-(sample,
chr) shuffle, so lowering the cap only removes reads. Background reads are a seeded sample of
`--background` raw reads whose IDs are not selected. gzip mtime is 0 for byte-identical rebuilds.

Usage:
  select_extract.py --sample S --bam PRE_DEDUP.bam --fastq RAW.fq[.gz] --chrs chrs.txt \
      --outdir DIR --seed 20261004 --cap 250 --background 150 [--umi-suffix 10]
"""
import argparse
import collections
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys


def norm_id(qname, umi_len):
    """Strip the `_<UMI>` suffix umi_tools extract appended to the read name."""
    if umi_len:
        head, sep, umi = qname.rpartition("_")
        if sep and len(umi) == umi_len:
            return head
    return qname


def rng(seed, *keys):
    """Independent, reproducible RNG per (seed, keys); not affected by PYTHONHASHSEED."""
    h = hashlib.sha256(":".join([str(seed)] + list(keys)).encode()).hexdigest()
    return random.Random(int(h[:16], 16))


def open_fq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def fq_records(path):
    with open_fq(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            rec = [h, fh.readline(), fh.readline(), fh.readline()]
            if not rec[3]:
                sys.exit("truncated FASTQ record in %s" % path)
            yield h[1:].split(None, 1)[0], "".join(rec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--chrs", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--cap", type=int, required=True)
    ap.add_argument("--background", type=int, required=True)
    ap.add_argument("--umi-suffix", type=int, default=0, help="UMI length appended to read names (0 = none)")
    a = ap.parse_args()

    chrs = [c.strip() for c in open(a.chrs) if c.strip()]
    keep = set(chrs)

    # --- §4.3 read-ID assignment from the pre-dedup BAM ---------------------------------
    primary_chr, first_sel = {}, {}
    p = subprocess.Popen(["samtools", "view", a.bam], stdout=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout:
        f = line.split("\t", 3)
        rid = norm_id(f[0], a.umi_suffix)
        flag = int(f[1])
        if not flag & 0x904:  # primary, mapped
            primary_chr[rid] = f[2]
        if f[2] in keep and rid not in first_sel:
            first_sel[rid] = f[2]
    if p.wait() != 0:
        sys.exit("samtools view failed on %s" % a.bam)
    by_chr = collections.defaultdict(list)
    for rid, c in first_sel.items():
        pc = primary_chr.get(rid)
        by_chr[pc if pc in keep else c].append(rid)

    capped = set()
    counts = {}
    for c in chrs:
        ids = sorted(by_chr.get(c, []))
        rng(a.seed, a.sample, c).shuffle(ids)
        capped.update(ids[: a.cap])
        counts[c] = {"selected": len(ids), "capped": min(len(ids), a.cap)}
    selected = set(first_sel)

    # --- §4.4 background: seeded reservoir sample of raw reads not selected --------------
    r = rng(a.seed, a.sample, "background")
    bg, n_other, n_raw = [], 0, 0
    for rid, _ in fq_records(a.fastq):
        n_raw += 1
        if rid in selected:
            continue
        n_other += 1
        if len(bg) < a.background:
            bg.append(rid)
        else:
            j = r.randrange(n_other)
            if j < a.background:
                bg[j] = rid
    bg = set(bg)

    # --- §4.5 extraction (original order, raw records) + check L2 ------------------------
    os.makedirs(os.path.join(a.outdir, "full"), exist_ok=True)
    os.makedirs(os.path.join(a.outdir, "capped"), exist_ok=True)
    os.makedirs(os.path.join(a.outdir, "ids"), exist_ok=True)
    want_full, want_cap = selected | bg, capped | bg
    out_full = gzip.GzipFile(os.path.join(a.outdir, "full", a.sample + ".fastq.gz"), "wb", 9, mtime=0)
    out_cap = gzip.GzipFile(os.path.join(a.outdir, "capped", a.sample + ".fastq.gz"), "wb", 9, mtime=0)
    found_full, found_cap, seen = 0, 0, set()
    for rid, rec in fq_records(a.fastq):
        if rid in seen:
            sys.exit("duplicate read ID in raw FASTQ: %s" % rid)
        if rid in want_full:
            seen.add(rid)
            out_full.write(rec.encode())
            found_full += 1
            if rid in want_cap:
                out_cap.write(rec.encode())
                found_cap += 1
    out_full.close()
    out_cap.close()

    for name, ids in (("selected", selected), ("capped", capped), ("background", bg)):
        with open(os.path.join(a.outdir, "ids", "%s.%s.txt" % (a.sample, name)), "w") as fh:
            fh.write("".join(i + "\n" for i in sorted(ids)))

    summary = {
        "sample": a.sample, "bam": a.bam, "fastq": a.fastq, "seed": a.seed, "cap": a.cap,
        "raw_reads": n_raw, "selected_reads": len(selected), "capped_reads": len(capped),
        "background_reads": len(bg), "per_chr": counts,
        "L2_full": {"wanted": len(want_full), "found": found_full},
        "L2_capped": {"wanted": len(want_cap), "found": found_cap},
    }
    with open(os.path.join(a.outdir, "ids", a.sample + ".summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    if found_full != len(want_full) or found_cap != len(want_cap):
        sys.exit("L2 FAILED for %s: %s" % (a.sample, json.dumps(summary)))
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
