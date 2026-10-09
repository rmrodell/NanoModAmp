#!/usr/bin/env python3
"""Golden plan §4.7 + §4.10: assemble tests/data/golden/ from the build directory.

Layout (matches conf/test_golden.config on wp/0-scaffold):
  tests/data/golden/
    README.md  provenance.json  MANIFEST.md5
    endogenous/      samplesheet_20250418.csv, samplesheet_20251022.csv, params_20250418.yaml,
                     params_20251022.yaml, params_merged.yaml, name_map.tsv, analyses.yaml (calling on
                     the counts-level merge of both runs, D32), fastq/, reference/, targets.bed, expected/
    mpra_invitro/    samplesheet.csv, params.yaml, analyses.yaml, fastq/, reference/, targets.bed, expected/
    mpra_incell/     same as mpra_invitro
  expected/ = paper_reference/ (reference only), legacy_full_selected/ (V1 counts),
              legacy_rerun/ (V2, binding), published_rerun/ (D30, reference)

Deterministic: FASTQs are byte copies of select/<run>/capped, references are bgzip'd (no
timestamps), generated text has no dates, MANIFEST.md5 is sorted. Rerunning gives an identical
MANIFEST.md5.

Usage: assemble_package.py [--golden DIR] [--out DIR] [--allow-missing]
"""
import argparse
import csv
import glob
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys

BUILD = "/scratch/users/rodell/NanoModAmp_build"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SEED, CAP, BACKGROUND, SIZE_CAP_MB = 20261004, 250, 150, 25.0

REF_ENDO = "/home/groups/nicolemm/rodell/fastas/endoG1G2.fasta"
REF_POOL1 = "/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta"
BED_ENDO = "/home/groups/nicolemm/rodell/fastas/endoG1G2_delpos.bed"
BED_POOL1 = "/home/groups/nicolemm/rodell/pool1/pool1_cleaned_delpos_noadapters.bed"
CHR_ENDO = BUILD + "/endoBID_data/endo_chr.txt"
CHR_MPRA = BUILD + "/pool1_chr.txt"
PAPER_COUNTS = {
    "endo_20250418": BUILD + "/endoBID_data/20250418_BIDdetect_data.txt",
    "endo_20251022": BUILD + "/endoBID_data/20251022_BIDdetect_data.txt",
    "mpra_invitro": "/oak/stanford/groups/nicolemm/rodell/BIDamplicon/Pool1/invitro/counts_delpos/BIDdetect_data_invitro_delpos.txt",
    "mpra_incell": "/oak/stanford/groups/nicolemm/rodell/BIDamplicon/Pool1/incell/counts_delpos/BIDdetect_data_incell_delpos.txt",
}
FIG3 = "/home/users/rodell/PUS7regulation2026/Figure3/plots"
PAPER_CALLS = {  # shipped paper site-call tables (PUS7regulation2026 Figure3/plots)
    "mpra_invitro": [FIG3 + "/invitro_modification_significance.tsv"],
    "mpra_incell": [FIG3 + "/WT_mod_Both_significance.tsv", FIG3 + "/PUS7_dep_union_significant_summary.tsv"],
}
VEC_ENDO = {"P36": "WT", "pLKO": "WT", "P97": "KD", "shPUS7": "KD", "WT": "WT", "KD": "KD"}

warnings = []


def warn(msg):
    warnings.append(msg)
    print("WARNING: " + msg, file=sys.stderr)


def need(path, allow_missing, what):
    if os.path.exists(path):
        return True
    if allow_missing:
        warn("missing %s: %s" % (what, path))
        return False
    sys.exit("ERROR: missing %s: %s (rerun with --allow-missing for a dry run)" % (what, path))


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        fh.write(text)


def chrs_of(path):
    return [l.strip() for l in open(path) if l.strip()]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# --- samples (approved G1-b sets, from the L1 inputs) -------------------------------------------
def samples(golden):
    l1 = os.path.join(golden, "l1", "inputs")
    out = {}
    out["endo_20250418"] = sorted(os.path.basename(l.strip())[:-3] for l in open(os.path.join(l1, "endo_20250418", "fastqs.txt")) if l.strip())
    for e in ("endo_20251022", "mpra_invitro"):
        out[e] = sorted(l.split(":")[0].strip() for l in open(os.path.join(l1, e, "sample_map.txt")) if l.strip())
    out["mpra_incell"] = sorted(os.path.basename(f)[:-3] for f in glob.glob(os.path.join(l1, "mpra_incell", "fq", "*.fq")))
    return out


def endo_row(run, s):
    """Original sample name -> sample-sheet row (vector WT/KD, run-unique rep as in published_rerun)."""
    ct, vec, rep, trt = s.split("_")
    new_rep = "%sr%s" % (run[4:], rep)
    sid = "%s_%s_%s_%s" % (ct, VEC_ENDO[vec], new_rep, trt)
    note = "original sample %s (run %s)" % (s, run)
    if vec != VEC_ENDO[vec]:
        note += "; vector %s -> %s" % (vec, VEC_ENDO[vec])
    return {"sample_id": sid, "treat": trt, "rep": new_rep, "celltype": ct, "vector": VEC_ENDO[vec], "run": run, "notes": note, "_orig": s}


def invitro_row(s):
    ct, vec, trt, rep = s.split("_")  # celltype_vector_treat_rep
    notes = []
    if trt == "in":
        trt = "input"
        notes.append("treat 'in' -> input")
    if vec == "noPUS" and trt == "BS":
        trt = "input"
        notes.append("noPUS BS labelled input (D14)")
    return {"sample_id": s, "treat": trt, "rep": rep, "celltype": ct, "vector": vec, "notes": "; ".join(notes), "_orig": s}


def incell_row(s):
    ct, vec, rep, trt = s.split("_")
    return {"sample_id": s, "treat": trt, "rep": rep, "celltype": ct, "vector": vec, "notes": "", "_orig": s}


def write_sheet(path, rows, cols):
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(cols)
    for r in sorted(rows, key=lambda r: r["sample_id"]):
        w.writerow([r[c] for c in cols])
    write(path, buf.getvalue())


# --- references and BED -----------------------------------------------------------------------
def add_reference(src, dst_dir, name):
    """bgzip (deterministic, no timestamp) + .fai/.gzi; fresh index (the pool1 .fai on disk is stale)."""
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, name + ".fa.gz")
    with open(dst, "wb") as out:
        subprocess.run(["bgzip", "-c", "-l", "9", src], stdout=out, check=True)
    for ext in (".fai", ".gzi"):
        if os.path.exists(dst + ext):
            os.remove(dst + ext)
    subprocess.run(["samtools", "faidx", dst], check=True)


def write_bed(src, chrs, dst):
    """Selected-target BED converted to standard 0-based (G1-e, bed_coordinates=bed0).

    The source delpos BEDs have start = end = the 1-based site; the package BED gets
    start = site - 1, end = site (all other columns unchanged). Fails if any selected row
    does not have start == end.
    """
    keep = set(chrs)
    rows = [l.rstrip("\r\n").split("\t") for l in open(src) if l.split("\t")[0] in keep]
    found = {r[0] for r in rows}
    if found != keep:
        sys.exit("ERROR: chrs missing from %s: %s" % (src, sorted(keep - found)))
    bad = [r[:3] for r in rows if r[1] != r[2]]
    if bad:
        sys.exit("ERROR: %s rows without start == end (expected single-site delpos BED): %s" % (src, bad))
    out = []
    for r in rows:
        site = int(r[2])
        out.append("\t".join([r[0], str(site - 1), str(site)] + r[3:]) + "\n")
    write(dst, "".join(out))
    return {(r[0], int(r[2])) for r in rows}  # 1-based sites


# --- analyses (treatment + PUS7 dependency only; cell-type-specific out of scope) ----------------
def yaml_analyses(analyses, site_sets, header):
    def val(v):
        if isinstance(v, list):
            return "[" + ", ".join(str(x) for x in v) + "]"
        if isinstance(v, dict):
            return "{" + ", ".join("%s: %s" % (k, val(x)) for k, x in v.items()) + "}"
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, str):
            return '"%s"' % v if (v == "" or v.startswith("#") or "|" in v or ":" in v) else v
        return str(v)
    out = [header, "analyses:"]
    for a in analyses:
        first = True
        for k, v in a.items():
            out.append(("  - " if first else "    ") + "%s: %s" % (k, val(v)))
            first = False
    if site_sets:
        out.append("site_sets:")
        for s in site_sets:
            first = True
            for k, v in s.items():
                out.append(("  - " if first else "    ") + "%s: %s" % (k, val(v)))
                first = False
    return "\n".join(out) + "\n"


def analyses_invitro():
    a = [{"name": "invitro", "type": "treatment", "subset": {"celltype": ["IV"]}, "random_effects": "",
          "sesoi": 0.05, "fdr": 0.05, "plot_all_sites": True, "colors": {"modified": "#c154c1", "input": "#eee8aa"}}]
    return yaml_analyses(a, [], "# Figure 3 in vitro (modification_analysis.R). noPUS samples are treat=input (D14).")


def analyses_cell(wt_vectors, factor_pairs, re_both_treat, re_single_treat, re_both_factor, header, union_name):
    a, union = [], []
    for cond in ("HepG2", "293T", "Both"):
        cts = ["HepG2", "293T"] if cond == "Both" else [cond]
        a.append({"name": "WT_mod_%s" % cond, "type": "treatment", "subset": {"celltype": cts, "vector": wt_vectors},
                  "random_effects": re_both_treat if cond == "Both" else re_single_treat, "sesoi": 0.05, "fdr": 0.05})
    for cond in ("HepG2", "293T", "Both"):
        cts = ["HepG2", "293T"] if cond == "Both" else [cond]
        for comp, levels in factor_pairs:
            name = "PUS7_dep_%s_%s" % (cond, comp)
            a.append({"name": name, "type": "factor", "subset": {"celltype": cts, "vector": levels}, "factor": "vector",
                      "levels": levels, "random_effects": re_both_factor if cond == "Both" else "", "direction": "positive",
                      "sesoi": 0.05, "fdr": 0.05})
            union.append(name)
    sets = [{"name": union_name, "op": "union", "of": union, "quartiles_by": {"analysis": "WT_mod_Both", "column": "delta_delrate"}}]
    return yaml_analyses(a, sets, header)


# --- expected outputs ---------------------------------------------------------------------------
def copy_tree(src, dst, allow_missing, what):
    if not need(src, allow_missing, what):
        return
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def paper_counts(path, meta_cols, keep_samples, sites, dst):
    """Paper count-table rows for included samples x selected sites (reference only)."""
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        hdr = next(r)
        idx = [hdr.index(c) for c in meta_cols]
        ci, pi = hdr.index("chr"), hdr.index("pos")
        rows = [row for row in r if tuple(row[i] for i in idx) in keep_samples and (row[ci], int(row[pi])) in sites]
    write(dst, "\t".join(hdr) + "\n" + "".join("\t".join(x) + "\n" for x in rows))
    return len(rows)


def paper_calls(path, chrs, dst):
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        hdr = next(r)
        ci = hdr.index("chr")
        rows = [row for row in r if row[ci] in chrs]
    write(dst, "\t".join(hdr) + "\n" + "".join("\t".join(x) + "\n" for x in rows))
    return len(rows)


# --- main ---------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--golden", default=BUILD + "/golden")
    ap.add_argument("--out", default=os.path.join(REPO, "tests", "data", "golden"))
    ap.add_argument("--allow-missing", action="store_true", help="dry run: warn instead of failing on missing inputs")
    a = ap.parse_args()
    G = a.golden
    if os.path.exists(a.out):
        shutil.rmtree(a.out)
    os.makedirs(a.out)

    S = samples(G)
    chr_endo, chr_mpra = chrs_of(CHR_ENDO), chrs_of(CHR_MPRA)
    chr_run = {"endo_20250418": [c for c in chr_endo if c.split("_")[0] in ("RHBDD2", "HDAC6")],
               "endo_20251022": [c for c in chr_endo if c.split("_")[0] == "STIM1"]}
    counts_summary = {}

    # ---------------- endogenous (two runs, merged at the counts level for calling, D32) --------
    E = os.path.join(a.out, "endogenous")
    rows_by_run, name_map = {}, ["run\toriginal_sample\tsample_id"]
    for run, key in (("20250418", "endo_20250418"), ("20251022", "endo_20251022")):
        rows = [endo_row(run, s) for s in S[key]]
        rows_by_run[run] = rows
        for r in rows:
            name_map.append("%s\t%s\t%s" % (run, r["_orig"], r["sample_id"]))
            src = os.path.join(G, "select", key, "capped", r["_orig"] + ".fastq.gz")
            r["fastq"] = "fastq/%s.fastq.gz" % r["sample_id"]
            if need(src, a.allow_missing, "capped FASTQ"):
                os.makedirs(os.path.join(E, "fastq"), exist_ok=True)
                shutil.copyfile(src, os.path.join(E, r["fastq"]))
        write_sheet(os.path.join(E, "samplesheet_%s.csv" % run), rows,
                    ["sample_id", "fastq", "treat", "rep", "celltype", "vector", "run", "notes"])
    write(os.path.join(E, "name_map.tsv"), "\n".join(name_map) + "\n")
    add_reference(REF_ENDO, os.path.join(E, "reference"), "endoG1G2")
    sites_endo = write_bed(BED_ENDO, chr_endo, os.path.join(E, "targets.bed"))
    common = "library_type: endogenous\nfasta: reference/endoG1G2.fa.gz\nbed: targets.bed\nbed_coordinates: bed0  # G1-e: targets.bed is standard 0-based\n"
    write(os.path.join(E, "params_20250418.yaml"),
          "# Run 20250418: reads lack the 5' ONT adapter -> opt-in 3'-only trim for THESE samples only (D31, R-29)\ninput: samplesheet_20250418.csv\n" + common + "ont_adapter_mode: three_prime_only\n")
    write(os.path.join(E, "params_20251022.yaml"), "# Run 20251022: pipeline defaults\ninput: samplesheet_20251022.csv\n" + common)
    write(os.path.join(E, "params_merged.yaml"),
          "# Site calling on the counts-level merge of both runs (D32): pass the counts_merged.tsv of the two runs above\n"
          "input_counts: <outdir_20250418>/counts/counts_merged.tsv,<outdir_20251022>/counts/counts_merged.tsv\nanalyses: analyses.yaml\n")
    write(os.path.join(E, "analyses.yaml"), analyses_cell(
        ["WT"], [("WT_v_KD", ["KD", "WT"])], "(1|celltype)", "", "(1|celltype)",
        "# Figure 2 endogenous under the Figure 3 method (D13), on the merged counts of runs 20250418 + 20251022 (D32).\n"
        "# Reps are run-unique (0418r1, 1022r2, ...); batch pairing holds within a run (D17).\n"
        "# Random effects (Becca, 2026-10-09, D33): replicate (1|rep), always in the model, plus (1|celltype) when both cell types are analysed together.", "PUS7_dep_union"))

    # ---------------- MPRA ----------------------------------------------------------------------
    mpra = {
        "mpra_invitro": ([invitro_row(s) for s in S["mpra_invitro"]], "# Figure 3 in vitro (May run): pipeline defaults\n", analyses_invitro()),
        "mpra_incell": ([incell_row(s) for s in S["mpra_incell"]],
                        "# Figure 3 in cellulo, run 20241114: library built without ONT adapters or UMIs -> opt-in modes for THESE samples only (D29, R-26).\n"
                        "# All other data, including new in-cellulo data, use the defaults (orientation_adapters=ont, umi=true).\n"
                        "orientation_adapters: pool\numi: false\n",
                        analyses_cell(["P102", "P4"], [("WT_v_KD", ["P101", "P102"]), ("WT_v_OE", ["P4", "P3"])],
                                      "(1|vector) + (1|celltype)", "(1|vector)", "(1|celltype)",
                                      "# Figure 3 in cellulo (incell_analysis.R): WT_mod + PUS7 dependency; cell-type-specific analyses out of scope.\n"
                                      "# Paper wt_vectors also included HepG2 'WT', which is not in the package (G1-b).", "PUS7_dep_union")),
    }
    sites_mpra = None
    for exp, (rows, extra, analyses) in mpra.items():
        D = os.path.join(a.out, exp)
        for r in rows:
            src = os.path.join(G, "select", exp, "capped", r["_orig"] + ".fastq.gz")
            r["fastq"] = "fastq/%s.fastq.gz" % r["sample_id"]
            if need(src, a.allow_missing, "capped FASTQ"):
                os.makedirs(os.path.join(D, "fastq"), exist_ok=True)
                shutil.copyfile(src, os.path.join(D, r["fastq"]))
        write_sheet(os.path.join(D, "samplesheet.csv"), rows, ["sample_id", "fastq", "treat", "rep", "celltype", "vector", "notes"])
        add_reference(REF_POOL1, os.path.join(D, "reference"), "pool1")
        sites_mpra = write_bed(BED_POOL1, chr_mpra, os.path.join(D, "targets.bed"))
        write(os.path.join(D, "params.yaml"), extra +
              "input: samplesheet.csv\nlibrary_type: mpra\nfasta: reference/pool1.fa.gz\nbed: targets.bed\nanalyses: analyses.yaml\n"
              "bed_coordinates: bed0  # G1-e: targets.bed is standard 0-based\n")
        write(os.path.join(D, "analyses.yaml"), analyses)

    # ---------------- expected outputs ----------------------------------------------------------
    exp_dirs = {"endogenous": "endogenous", "mpra_invitro": "mpra_invitro", "mpra_incell": "mpra_incell"}
    for exp in exp_dirs:
        X = os.path.join(a.out, exp, "expected")
        for kind in ("legacy_full_selected", "legacy_rerun", "published_rerun"):
            copy_tree(os.path.join(G, "expected", kind, exp), os.path.join(X, kind), a.allow_missing, "expected/%s/%s" % (kind, exp))
    # paper_reference: counts (included samples x selected sites) + shipped site calls
    keep = {
        "endo_20250418": {tuple(s.split("_")) for s in S["endo_20250418"]},
        "endo_20251022": {tuple(s.split("_")) for s in S["endo_20251022"]},
        "mpra_invitro": {tuple(s.split("_")) for s in S["mpra_invitro"]},
        "mpra_incell": {tuple(s.split("_")) for s in S["mpra_incell"]},
    }
    meta = {"endo_20250418": ["celltype", "vector", "rep", "treat"], "endo_20251022": ["celltype", "vector", "rep", "treat"],
            "mpra_invitro": ["celltype", "vector", "treat", "rep"], "mpra_incell": ["celltype", "vector", "rep", "treat"]}
    n_ref = {}
    for key, path in PAPER_COUNTS.items():
        exp = "endogenous" if key.startswith("endo") else key
        sites = sites_endo if exp == "endogenous" else sites_mpra
        dst = os.path.join(a.out, exp, "expected", "paper_reference",
                           ("counts_%s.tsv" % key.split("_")[1]) if exp == "endogenous" else "counts.tsv")
        if need(path, a.allow_missing, "paper counts"):
            n_ref[key] = paper_counts(path, meta[key], keep[key], sites, dst)
    for exp, paths in PAPER_CALLS.items():
        for p in paths:
            if need(p, a.allow_missing, "paper site calls"):
                n_ref["%s:%s" % (exp, os.path.basename(p))] = paper_calls(p, set(chr_mpra), os.path.join(a.out, exp, "expected", "paper_reference", os.path.basename(p)))

    # ---------------- provenance.json -----------------------------------------------------------
    inputs = []
    for tsv in sorted(glob.glob(os.path.join(G, "provenance", "*_inputs.tsv"))):
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            inputs.append({"file": os.path.basename(r["path"]), "dir": os.path.basename(os.path.dirname(r["path"])),
                           "bytes": int(r["bytes"]), "md5": r.get("md5") or None})
    inputs.sort(key=lambda d: (d["dir"], d["file"]))
    per_sample = []
    summ = os.path.join(G, "select", "select_summary.tsv")
    if need(summ, a.allow_missing, "selection summary"):
        per_sample = sorted(csv.DictReader(open(summ), delimiter="\t"), key=lambda d: (d["experiment"], d["sample"], d["chr"]))
    versions = {}
    for name in ("legacy", "published"):
        p = os.path.join(REPO, "dev", name, "versions.txt")
        if need(p, a.allow_missing, "versions.txt"):
            versions[name] = [l.rstrip("\n") for l in open(p) if l.strip()]
    git = lambda repo: subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    prov = {
        "package": "NanoModAmp golden test package (golden_test_package_plan.md)",
        "paper_repo": {"PUS7regulation2026": git("/home/users/rodell/PUS7regulation2026"), "BIDamplicon": git("/home/users/rodell/BIDamplicon")},
        "build": {"seed": SEED, "cap_selected_reads_per_sample_target": CAP, "background_reads_per_sample": BACKGROUND,
                  "size_cap_mb": SIZE_CAP_MB, "umi_dedup_seed_legacy": SEED},
        "references": {"endogenous": os.path.basename(REF_ENDO), "mpra": os.path.basename(REF_POOL1)},
        "targets": {"endogenous": chr_endo, "mpra": chr_mpra},
        "samples": {k: v for k, v in S.items()},
        "per_sample_per_target_reads": per_sample,
        "paper_reference_rows": n_ref,
        "tool_versions": versions,
        "inputs": inputs,
        "warnings": warnings,
    }
    write(os.path.join(a.out, "provenance.json"), json.dumps(prov, indent=1, sort_keys=True) + "\n")

    # ---------------- README --------------------------------------------------------------------
    write(os.path.join(a.out, "README.md"), README)

    # ---------------- size check + manifest -----------------------------------------------------
    files = sorted(os.path.relpath(os.path.join(d, f), a.out) for d, _, fs in os.walk(a.out) for f in fs if f != "MANIFEST.md5")
    sizes = {}
    for f in files:
        top = f.split(os.sep)[0] if os.sep in f else "(root)"
        sizes[top] = sizes.get(top, 0) + os.path.getsize(os.path.join(a.out, f))
    total = sum(sizes.values()) / 1e6
    report = "\n".join("%-14s %8.2f MB" % (k, v / 1e6) for k, v in sorted(sizes.items())) + "\n%-14s %8.2f MB (cap %.0f MB)\n" % ("TOTAL", total, SIZE_CAP_MB)
    print(report, end="")
    write(os.path.join(a.out, "MANIFEST.md5"), "".join("%s  %s\n" % (md5(os.path.join(a.out, f)), f) for f in files))
    if warnings:
        print("%d warning(s); package is INCOMPLETE" % len(warnings), file=sys.stderr)
    if total > SIZE_CAP_MB:
        sys.exit("ERROR: package %.2f MB exceeds the %.0f MB cap" % (total, SIZE_CAP_MB))


README = """# NanoModAmp golden test package

A small public subset of the paper data (*RNA sequence, structure, and cell type specific
features drive pseudouridylation by PUS7*, Figures 2 and 3). Run it to check that an
installation reproduces known results. CI uses it as a regression test.

Built by `dev/golden/assemble_package.py` (plan: `docs/plans/golden_test_package_plan.md`;
decisions: `docs/DECISIONS.md`). `provenance.json` records inputs (basenames + md5), tool
versions, seed and per-sample read counts; `MANIFEST.md5` lists every file.

## Contents

| Directory | Data | Samples | Targets |
|---|---|---|---|
| `endogenous/` | Figure 2, runs 20250418 and 20251022 | 32 (16 per run) | RHBDD2, HDAC6 (20250418); STIM1 (20251022) |
| `mpra_invitro/` | Figure 3 in vitro, May run | 4 | 6 pool1 oligos (`pool1_chr.txt`) |
| `mpra_incell/` | Figure 3 in cellulo, run 20241114 | 32 | same 6 oligos |

Each directory has `fastq/` (raw reads, adapters intact, original order), `reference/` (the
**full** reference, bgzip + index), `targets.bed` (paper single-site BED filtered to the targets,
**converted to standard 0-based**: start = site − 1, end = site; run with `bed_coordinates=bed0`), sample sheet(s), `params*.yaml`, `analyses.yaml` and `expected/`.

## How the reads were chosen (G1-b, G1-d)

All reads with any alignment on a target in the legacy rerun's pre-dedup BAM, with their PCR
duplicates, then at most **250 reads per sample and target** (seeded, nested) plus **150 random
background reads per sample** (seeded) so the adapter, length and mapping filters are exercised.
Targets below the cap are kept whole. Seed 20261004.

## Running it

    nextflow run . -profile test_golden,<docker|singularity> --golden_experiment <name> --outdir <dir>

`<name>` is `endogenous_20250418`, `endogenous_20251022`, `mpra_invitro` or `mpra_incell`.
Endogenous site calling runs on the counts-level merge of both runs (D32):
`--input_counts <dir_20250418>/counts/counts_merged.tsv,<dir_20251022>/counts/counts_merged.tsv
--analyses tests/data/golden/endogenous/analyses.yaml` (see `endogenous/params_merged.yaml`).

### Non-default modes — these samples only

Pipeline defaults (ONT trim with both adapters, UMI extraction, deduplication) apply to all
data, including in vitro and in cellulo. Two paper libraries need opt-in settings, set only for
them in `conf/test_golden.config` and the `params*.yaml` here:

- **Run 20250418 (endogenous):** reads lack the 5′ ONT adapter → `ont_adapter_mode=three_prime_only` (D31, R-29).
- **Run 20241114 (in cellulo):** library built without ONT adapters or UMIs → `orientation_adapters=pool`, `umi=false` (D29, R-26).

### Sample-sheet conventions

- Endogenous vectors are normalised: P36/pLKO → WT, P97/shPUS7 → KD; reps are run-unique
  (`0418r1`, `1022r2`, …) because batch pairing holds only within a run. `endogenous/name_map.tsv`
  maps package names to the original sample names.
- In vitro noPUS samples are `treat=input` (D14), noted in `notes`.

## Expected outputs: what is binding

| Directory | Source | Status |
|---|---|---|
| `expected/legacy_rerun/` | Scripts that actually produced the paper numbers (`dev/legacy/`), rerun on this package (V2) | **Binding for counts.** Endogenous calls are reference only (legacy used `glmer`; the pipeline uses `bglmer`, D13). |
| `expected/published_rerun/` | Published Figure 2/3 scripts (`dev/published/`), rerun on this package (D30) | Reference; isolates implementation differences in WP7 |
| `expected/legacy_full_selected/` | Legacy counts on the selection before the read cap (V1) | Reference; V1 matched the paper counts |
| `expected/paper_reference/` | Paper count tables and shipped site-call tables, restricted to these samples and targets | Reference only: paper calls used all samples and sites |

Notes:
- The paper processing differs from the published scripts (R-25): minimap2 `-a -k5`, no
  primary/MAPQ filter, `umi_tools dedup --method directional`. The legacy harness adds a fixed
  `--random-seed` (G1-c), so `legacy_rerun` is deterministic but its dedup representatives
  differ from the paper's.
- Paper counts for 20250418 were made on **pre-dedup** BAMs; 20251022 and in vitro on dedup BAMs;
  in cellulo has no dedup (R-27). V1 was exact for non-deduplicated counts and within tolerance
  for deduplicated ones (G1-c).
- `endoG1G2.fasta` is used for all endogenous runs (F2); target sequences are identical to the
  per-run references the paper used.
- BED convention (R-15, gate G1-e): the pipeline default is `bed_coordinates=bed0`. The paper's
  delpos BEDs have start = end = the 1-based site (the legacy counter adds no +1); the package
  `targets.bed` is converted to standard 0-based. The legacy and published reruns used the
  original start = end BED with the legacy code. Positions in all expected outputs are 1-based
  sites and are identical under both conventions.
- Run 20250418 paper counts also contain window-start rows from its whole-amplicon BED; these
  are a documented paper artifact (`dev/golden/R15_bed_convention.md`) and are not in the package.
- Known gaps (G1-b): no Inconclusive, input-delrate = 0 or Modified-not-PUS7-dependent MPRA
  site; cell-type-specific analyses are out of scope.
"""

if __name__ == "__main__":
    main()
