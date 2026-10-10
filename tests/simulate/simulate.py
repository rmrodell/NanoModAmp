#!/usr/bin/env python3
"""Synthetic Nano-BID-Amp data and truth tables (plan §7 WP1).

Generates, for every experiment in the design YAML: per-sample FASTQs, a sample sheet (§5.1),
params and analyses YAML (§5.4), and truth tables:

  truth/funnel.tsv          expected records per read-funnel step (§5.2) and sample
  truth/both_orientations.tsv  n_reads_in_both_orientations per sample
  truth/site_counts.tsv     expected per-site counts (§5.3 columns) for the single-site BED
  truth/amplicon_counts.tsv expected counts at every T of every target (whole-target BED, D11)
  truth/categories.tsv      expected calls for clear-cut sites, per analysis
  truth/reads.tsv.gz        every read: class, molecule, orientation and fate

References and BEDs (standard 0-based, `bed_coordinates=bed0`) go to <outdir>/reference/.

Truth for trimming is computed by emulating the pipeline's cutadapt passes on the actual read
sequences (adapters are inserted without errors, so exact matching is what cutadapt does).
Mapping outcomes are by construction: every trimmed record either equals a known (edited) target
sequence, which maps uniquely with MAPQ >= 30, comes from the near-duplicate pair (MAPQ 0), or is
random sequence (unmapped). tests/simulate/test_tools.py checks these assumptions with minimap2.
Deletions sit on a T flanked by non-T bases and inserted bases differ from both neighbours, so
every indel has exactly one alignment.

Deterministic: the same design and seed give byte-identical output (gzip mtime 0, sorted keys).

Usage: simulate.py --design design.yaml --outdir tests/data/synthetic
"""
import argparse
import csv
import gzip
import hashlib
import io
import os
import random
import shutil
import sys

import yaml

BASES = "ACGT"
COMP = str.maketrans("ACGTN", "TGCAN")
COUNT_COLS = ["chr", "pos", "gene", "totalReads", "A.count", "C.count", "G.count", "T.count",
              "Deletion.count", "Insertion.count", "ref", "kmer", "strand", "delrate"]


def rc(s):
    return s.translate(COMP)[::-1]


def rng_for(seed, *keys):
    """Independent RNG per purpose; stable across Python versions and PYTHONHASHSEED."""
    h = hashlib.sha256(":".join(str(k) for k in (seed,) + keys).encode()).hexdigest()
    return random.Random(int(h[:16], 16))


def rround(x):
    return int(x + 0.5)


# --------------------------------------------------------------------------------------------
# References
# --------------------------------------------------------------------------------------------
def random_seq(r, n, avoid=None):
    """Random sequence without homopolymers longer than 3."""
    out = []
    while len(out) < n:
        b = r.choice(BASES)
        if len(out) >= 3 and out[-1] == out[-2] == out[-3] == b:
            continue
        out.append(b)
    s = "".join(out)
    if avoid:
        for a in avoid:
            if a in s or rc(a) in s:
                return random_seq(r, n, avoid)
    return s


def build_reference(design, name):
    spec = design["references"][name]
    r = rng_for(design["seed"], "reference", name, spec.get("seed_offset", 0))
    adapters = list(design["adapters"].values())
    targets, seqs = [], {}
    for t in spec["targets"]:
        if t.get("near_duplicate_of"):
            base = seqs[t["near_duplicate_of"]]
            seq = base + random_seq(r, t["length"] - len(base))
        else:
            seq = list(random_seq(r, t["length"], adapters))
            for s in t["sites"]:
                i = s["pos"] - 1
                if not (20 <= s["pos"] <= t["length"] - 20):
                    sys.exit("site %s:%d must be >= 20 nt from both ends" % (t["name"], s["pos"]))
                seq[i] = "T"
                # Flanks are never T, so a single-T deletion has one alignment
                seq[i - 1] = r.choice("ACG")
                seq[i + 1] = r.choice("ACG")
            seq = "".join(seq)
        sites = []
        for s in t["sites"]:
            ins_base = None
            if design["site_types"][s["type"]].get("p_ins"):
                right = seq[s["pos"]]          # base after the T (0-based index pos)
                ins_base = sorted(set("ACG") - {right})[0]
            sites.append({"pos": s["pos"], "type": s["type"], "ins_base": ins_base})
        seqs[t["name"]] = seq
        targets.append({"name": t["name"], "gene": t["gene"], "seq": seq, "sites": sites})
    return targets


def write_reference(outdir, name, targets):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, name + ".fa"), "w") as fh:
        for t in targets:
            fh.write(">%s\n" % t["name"])
            for i in range(0, len(t["seq"]), 60):
                fh.write(t["seq"][i:i + 60] + "\n")
    with open(os.path.join(outdir, name + ".sites.bed"), "w") as fh:
        for t in targets:
            for s in t["sites"]:
                fh.write("%s\t%d\t%d\t%s\t0\t+\n" % (t["name"], s["pos"] - 1, s["pos"], t["gene"]))
    with open(os.path.join(outdir, name + ".amplicons.bed"), "w") as fh:
        for t in targets:
            fh.write("%s\t0\t%d\t%s\t0\t+\n" % (t["name"], len(t["seq"]), t["gene"]))


# --------------------------------------------------------------------------------------------
# Molecules and reads
# --------------------------------------------------------------------------------------------
def site_p(spec, key, treat, vector, rep, mult):
    v = (spec.get(key) or {}).get(treat, 0.0)
    if isinstance(v, dict):
        v = v[vector]
    if isinstance(v, list):
        return v[rep - 1]
    return v * mult


def edited_insert(seq, outcomes, sites):
    """Apply per-site outcomes (T, del, ins, C) to a target sequence, right to left."""
    s = list(seq)
    for site in sorted(sites, key=lambda x: -x["pos"]):
        o = outcomes[site["pos"]]
        i = site["pos"] - 1
        if o == "del":
            del s[i]
        elif o == "ins":
            s.insert(i + 1, site["ins_base"])
        elif o == "C":
            s[i] = "C"
    return "".join(s)


def make_umis(r, n, length, min_dist):
    umis = []
    tries = 0
    while len(umis) < n:
        u = random_seq(r, length)
        tries += 1
        if tries > 200000:
            sys.exit("could not draw %d UMIs with Hamming distance >= %d" % (n, min_dist))
        if all(sum(a != b for a, b in zip(u, v)) >= min_dist for v in umis):
            umis.append(u)
    return umis


def one_mismatch(r, umi):
    i = r.randrange(len(umi))
    return umi[:i] + r.choice([b for b in BASES if b != umi[i]]) + umi[i + 1:]


class Builder:
    """Read layout per library mode (adapter and UMI placement)."""

    def __init__(self, design, params):
        a = design["adapters"]
        self.S5, self.S3 = a["ont_sense_5p"], a["ont_sense_3p"]
        self.P5, self.P3 = a["pool_5p"], a["pool_3p"]
        self.mpra = params["library_type"] == "mpra"
        self.pool = params["orientation_adapters"] == "pool"
        self.three = params["ont_adapter_mode"] == "three_prime_only"
        self.umi = params["umi"]

    def sense(self, insert, umi):
        """Sense-orientation molecule as sequenced."""
        if self.pool:
            return self.P5 + insert + self.P3
        core = (self.P5 + insert + self.P3) if self.mpra else insert
        core += umi if self.umi else ""
        if self.three:
            return core + self.S3
        return self.S5 + core + self.S3


def build_sample(design, exp, targets, sample, builder):
    seed = design["seed"]
    r = rng_for(seed, exp["name"], sample["sample_id"], "molecules")
    umi_len, min_d = design["umi_length"], design["min_umi_hamming"]
    mult = exp["rep_multiplier"][sample["rep"]]
    copies_w = sorted(exp["copies"].items())
    molecules, reads = [], []
    for t in targets:
        n = exp["molecules_override"].get(t["name"], exp["molecules_per_target"])
        if isinstance(n, dict):
            n = n[sample["rep"]]
        if n == 0:
            continue
        outcomes = [dict() for _ in range(n)]
        for s in t["sites"]:
            spec = design["site_types"][s["type"]]
            nd = rround(site_p(spec, "p_del", sample["treat"], sample["vector"], sample["rep"], mult) * n)
            ni = rround(site_p(spec, "p_ins", sample["treat"], sample["vector"], sample["rep"], mult) * n)
            ns = rround(site_p(spec, "p_sub", sample["treat"], sample["vector"], sample["rep"], mult) * n)
            if nd + ni + ns > n:
                sys.exit("site %s:%d over-assigned" % (t["name"], s["pos"]))
            idx = list(range(n))
            rng_for(seed, exp["name"], sample["sample_id"], t["name"], s["pos"]).shuffle(idx)
            labels = ["del"] * nd + ["ins"] * ni + ["C"] * ns + ["T"] * (n - nd - ni - ns)
            for i, lab in zip(idx, labels):
                outcomes[i][s["pos"]] = lab
        umis = make_umis(r, n, umi_len, min_d)
        for k in range(n):
            mid = "%s:%d" % (t["name"], k)
            ncopy = r.choices([c for c, _ in copies_w], weights=[w for _, w in copies_w])[0]
            mismatch = ncopy >= 3 and r.random() < exp["umi_mismatch_fraction"]
            insert = edited_insert(t["seq"], outcomes[k], t["sites"])
            mol = {"id": mid, "target": t["name"], "umi": umis[k], "outcomes": outcomes[k],
                   "insert": insert, "copies": ncopy, "multimapper": t["name"].startswith("NEARDUP")}
            molecules.append(mol)
            for c in range(ncopy):
                u = one_mismatch(r, umis[k]) if (mismatch and c == ncopy - 1) else umis[k]
                reads.append({"class": "multimapper" if mol["multimapper"] else "valid",
                              "parts": [(mol, u)], "umi_variant": u != umis[k]})
    # junk
    j = exp["junk"]
    jr = rng_for(seed, exp["name"], sample["sample_id"], "junk")
    valid_mols = [m for m in molecules if not m["multimapper"]]
    first = targets[0]
    for cls in sorted(j):
        for _ in range(j[cls]):
            umi = random_seq(jr, umi_len)
            if cls == "chimera":
                m1, m2 = jr.sample(valid_mols, 2)
                reads.append({"class": "chimera", "parts": [(m1, m1["umi"]), (m2, m2["umi"])],
                              "umi_variant": False})
                continue
            reads.append({"class": cls, "seq": junk_seq(cls, jr, builder, first, umi),
                          "parts": [], "umi_variant": False})
    return molecules, reads


def junk_seq(cls, r, b, target, umi):
    """Junk read classes (plan WP1); all built from the target set's first sequence or random."""
    ins = target["seq"]
    if cls == "no_adapter":
        return random_seq(r, 200)
    if cls == "adapter_dimer":
        if b.pool:
            return b.P5 + b.P3
        return (umi + b.S3) if b.three else (b.S5 + umi + b.S3)
    if cls == "unmappable":
        rnd = random_seq(r, 200 if not b.mpra else 130)
        return b.sense(rnd, umi)
    if cls == "too_short":
        frag = ins[:70] if (b.mpra and not b.pool) else ins[:100] if b.pool else ins[:80]
        return b.sense(frag, umi)
    if cls == "short_pool":                     # passes ONT trim, fails the pool trim length
        return b.sense(ins[:105], umi)
    if cls == "one_adapter_5p":
        full = b.sense(ins, umi)
        return full[:len(full) - len(b.P3 if b.pool else b.S3)]
    if cls == "one_adapter_3p":
        full = b.sense(ins, umi)
        return full[len(b.P5 if b.pool else b.S5):]
    sys.exit("unknown junk class " + cls)


def assemble_reads(design, exp, sample, builder, reads):
    """Turn read specs into sequences with orientation; assign IDs; shuffle order."""
    r = rng_for(design["seed"], exp["name"], sample["sample_id"], "orient")
    for rd in reads:
        if rd["parts"] and rd["class"] != "chimera":
            mol, u = rd["parts"][0]
            s = builder.sense(mol["insert"], u)
            if builder.three:
                rd["orientation"] = "sense"     # D31 libraries: sense reads only (documented)
            else:
                rd["orientation"] = r.choice(["sense", "antisense"])
            rd["seq"] = s if rd["orientation"] == "sense" else rc(s)
        elif rd["class"] == "chimera":
            (m1, u1), (m2, u2) = rd["parts"]
            rd["seq"] = builder.sense(m1["insert"], u1) + rc(builder.sense(m2["insert"], u2))
            rd["orientation"] = "both"
        else:
            rd["orientation"] = "sense"
            if r.random() < 0.5 and rd["class"] != "no_adapter" and not builder.three:
                rd["seq"] = rc(rd["seq"])
                rd["orientation"] = "antisense"
    order = list(range(len(reads)))
    r.shuffle(order)
    tag = hashlib.sha256(("%s:%s" % (exp["name"], sample["sample_id"])).encode()).hexdigest()[:8]
    out = []
    for n, i in enumerate(order):
        rd = reads[i]
        rd["id"] = "%s%06x" % (tag, n)
        out.append(rd)
    return out


# --------------------------------------------------------------------------------------------
# Pipeline emulation for truth
# --------------------------------------------------------------------------------------------
def linked_trim(seq, a5, a3, min_len):
    i = seq.find(a5)
    if i < 0:
        return None
    j = seq.find(a3, i + len(a5))
    if j < 0:
        return None
    out = seq[i + len(a5):j]
    return out if len(out) >= min_len else None


def emulate(design, builder, rd, registry):
    """Records per funnel step for one read, and the fate of each surviving record.

    Each trimmed record is traced back to the molecule it came from (sense output = first part,
    antisense output = last part) and checked against that molecule's insert and UMI.
    """
    p = design["pipeline"]
    steps = {}
    recs = []                                      # (sequence, part index)
    if builder.three:
        j = rd["seq"].find(builder.S3)
        out = rd["seq"][:j] if j >= 0 else rd["seq"]
        recs = [(out, 0)] if len(out) >= p["trim1_min_length"] else []
        steps["trim1_3prime"] = len(recs)
    else:
        a5, a3 = (builder.P5, builder.P3) if builder.pool else (builder.S5, builder.S3)
        sense = linked_trim(rd["seq"], a5, a3, p["trim1_min_length"])
        anti = linked_trim(rd["seq"], rc(a3), rc(a5), p["trim1_min_length"])
        steps["trim1_sense"] = int(sense is not None)
        steps["trim1_antisense"] = int(anti is not None)
        steps["antisense_rc"] = steps["trim1_antisense"]
        if sense is not None:
            recs.append((sense, 0))
        if anti is not None:
            recs.append((rc(anti), len(rd["parts"]) - 1))
        steps["merged"] = len(recs)
        rd["both_orientations"] = int(sense is not None and anti is not None)
    umis = [None] * len(recs)
    if builder.umi:
        umis = [x[-design["umi_length"]:] for x, _ in recs]
        recs = [(x[:-design["umi_length"]], k) for x, k in recs]
        steps["umi_extracted"] = len(recs)
    if builder.mpra and not builder.pool:
        kept = []
        for (x, k), u in zip(recs, umis):
            y = linked_trim(x, builder.P5, builder.P3, p["trim2_min_length"])
            if y is not None:
                kept.append(((y, k), u))
        recs, umis = [r_[0] for r_ in kept], [r_[1] for r_ in kept]
        steps["trim2_pool"] = len(recs)
    steps["primary"] = len(recs)
    fates = []
    for (x, k), u in zip(recs, umis):
        if not rd["parts"]:
            if x in registry or rd["class"] not in ("unmappable", "no_adapter"):
                raise RuntimeError("junk read %s (%s) reached mapping as a target sequence" % (rd["id"], rd["class"]))
            fates.append(("unmapped", None, u, x))
            continue
        mol, umi = rd["parts"][k]
        if x != mol["insert"] or (builder.umi and u != umi):
            raise RuntimeError("read %s: trimmed record does not match molecule %s" % (rd["id"], mol["id"]))
        fates.append(("mapq0" if mol["multimapper"] else "mapped", mol, u, x))
    steps["mapq_filtered"] = sum(f[0] == "mapped" for f in fates)
    return steps, fates


FUNNEL_ORDER = ["raw", "trim1_sense", "trim1_antisense", "antisense_rc", "merged", "trim1_3prime",
                "umi_extracted", "trim2_pool", "mapped", "primary", "mapq_filtered", "dedup"]


# --------------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------------
def gz_write(path, text):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=raw, mtime=0) as g:
            g.write(text.encode())


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def fmt(x):
    return repr(float(x)) if isinstance(x, float) else str(x)


def count_rows(targets, records_by_target, site_only):
    """Expected count rows from counted units (molecules or read records) per target."""
    rows = []
    for t in targets:
        units = records_by_target.get(t["name"], [])
        if not units:
            continue
        sites = {s["pos"]: s for s in t["sites"]}
        positions = sorted(sites) if site_only else [i + 1 for i, b in enumerate(t["seq"]) if b == "T"]
        for pos in positions:
            c = {"A": 0, "C": 0, "G": 0, "T": 0, "Deletion": 0, "Insertion": 0}
            for mol in units:
                o = mol["outcomes"].get(pos, "T")
                if o == "del":
                    c["Deletion"] += 1
                elif o == "ins":
                    c["T"] += 1
                    c["Insertion"] += 1
                else:
                    c[o] += 1
            total = sum(c.values())
            seq = t["seq"]
            kmer = seq[pos - 3:pos + 2] if 3 <= pos <= len(seq) - 2 else "NA"
            rows.append([t["name"], pos, t["gene"], total, c["A"], c["C"], c["G"], c["T"],
                         c["Deletion"], c["Insertion"], "T", kmer, "+", c["Deletion"] / total])
    return rows


def experiment_samples(exp):
    sm = exp["samples"]
    return [{"sample_id": "SIM_%s_%d_%s" % (v, rep, tr), "vector": v, "treat": tr, "rep": rep,
             "celltype": sm["celltype"]}
            for v in sm["vector"] for rep in sm["rep"] for tr in sm["treat"]]


def mapping_inputs(design, ename, sample_index=0):
    """Trimmed records that reach the aligner for one sample, with their expected fate.

    Used by test_tools.py to check the mapping assumptions behind the truth tables.
    Returns (targets, params, [(record_id, sequence, fate, molecule or None)]).
    """
    exp = dict(design["experiments"][ename], name=ename)
    targets = build_reference(design, exp["reference"])
    builder = Builder(design, exp["params"])
    s = experiment_samples(exp)[sample_index]
    molecules, reads = build_sample(design, exp, targets, s, builder)
    reads = assemble_reads(design, exp, s, builder, reads)
    registry = {m["insert"]: m for m in molecules}
    out = []
    for rd in reads:
        _, fates = emulate(design, builder, rd, registry)
        for k, (f, mol, _, seq) in enumerate(fates):
            out.append(("%s.%d" % (rd["id"], k), seq, f, mol))
    return targets, exp["params"], out


def simulate(design, outdir, params_root="tests/data/synthetic"):
    if design.get("error_model", "none") != "none":
        sys.exit("error_model %r is not implemented; use 'none' (deterministic tests)" % design["error_model"])
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    refdir = os.path.join(outdir, "reference")
    refs = {}
    for name in sorted(design["references"]):
        refs[name] = build_reference(design, name)
        write_reference(refdir, name, refs[name])

    for ename in sorted(design["experiments"]):
        exp = dict(design["experiments"][ename], name=ename)
        targets = refs[exp["reference"]]
        params = exp["params"]
        builder = Builder(design, params)
        edir = os.path.join(outdir, ename)
        os.makedirs(os.path.join(edir, "fastq"))
        os.makedirs(os.path.join(edir, "truth"))
        samples = experiment_samples(exp)
        funnel, both, site_rows, amp_rows, read_rows, prededup_rows = [], [], [], [], [], []
        for s in samples:
            molecules, reads = build_sample(design, exp, targets, s, builder)
            reads = assemble_reads(design, exp, s, builder, reads)
            registry = {}
            for m in molecules:
                if m["insert"] in registry and registry[m["insert"]]["target"] != m["target"]:
                    raise RuntimeError("two targets share an edited insert")
                registry.setdefault(m["insert"], m)
            # Near-duplicate reads carry NEARDUP_A's sequence; registry maps it to the multimapper
            totals = {k: 0 for k in FUNNEL_ORDER}
            totals["raw"] = len(reads)
            counted_reads, dedup_mols = {}, {}
            nboth = 0
            for rd in reads:
                steps, fates = emulate(design, builder, rd, registry)
                for k, v in steps.items():
                    totals[k] += v
                nboth += rd.get("both_orientations", 0)
                for f, mol, u, _ in fates:
                    if f == "mapped":
                        counted_reads.setdefault(mol["target"], []).append(mol)
                        dedup_mols.setdefault(mol["target"], {})[mol["id"]] = mol
                read_rows.append([s["sample_id"], rd["id"], rd["class"], rd["orientation"],
                                  ";".join(m["id"] for m, _ in rd["parts"]) or "-",
                                  ";".join(u for _, u in rd["parts"]) or "-",
                                  int(rd["umi_variant"]),
                                  ";".join(f[0] for f in fates) or "dropped_trim"])
            totals["dedup"] = sum(len(v) for v in dedup_mols.values())
            steps_present = ["raw"]
            if builder.three:
                steps_present.append("trim1_3prime")
            else:
                steps_present += ["trim1_sense", "trim1_antisense", "antisense_rc", "merged"]
            if builder.umi:
                steps_present.append("umi_extracted")
            if builder.mpra and not builder.pool:
                steps_present.append("trim2_pool")
            steps_present += ["mapped", "primary", "mapq_filtered"]
            if builder.umi:
                steps_present.append("dedup")
            for st in steps_present:
                # 'mapped' = records written by minimap2: one primary-or-unmapped record per read plus
                # secondary alignments, which depend on the aligner; truth gives the read count.
                val = totals["primary"] if st == "mapped" else totals[st]
                funnel.append([s["sample_id"], st, val])
            if not builder.three:
                both.append([s["sample_id"], nboth])
            units = {k: list(v.values()) for k, v in dedup_mols.items()} if builder.umi else counted_reads
            level = "dedup" if builder.umi else "reads"
            for row in count_rows(targets, units, True):
                site_rows.append([s["sample_id"], level] + row + [int(row[3] > design["pipeline"]["min_coverage"])])
            for row in count_rows(targets, units, False):
                amp_rows.append([s["sample_id"], level] + row + [int(row[3] > design["pipeline"]["min_coverage"])])
            for row in count_rows(targets, counted_reads, True):
                prededup_rows.append([s["sample_id"], "reads"] + row)
            fq = io.StringIO()
            q = design["quality_char"]
            for rd in reads:
                fq.write("@%s\n%s\n+\n%s\n" % (rd["id"], rd["seq"], q * len(rd["seq"])))
            gz_write(os.path.join(edir, "fastq", s["sample_id"] + ".fastq.gz"), fq.getvalue())
            s["fastq"] = "fastq/%s.fastq.gz" % s["sample_id"]

        t = os.path.join(edir, "truth")
        write_tsv(os.path.join(t, "funnel.tsv"), ["sample_id", "step", "records"], funnel)
        if both:
            write_tsv(os.path.join(t, "both_orientations.tsv"), ["sample_id", "n_reads_in_both_orientations"], both)
        hdr = ["sample_id", "level"] + COUNT_COLS + ["passes_min_coverage"]
        write_tsv(os.path.join(t, "site_counts.tsv"), hdr, [[fmt(x) for x in r] for r in site_rows])
        write_tsv(os.path.join(t, "amplicon_counts.tsv"), hdr, [[fmt(x) for x in r] for r in amp_rows])
        if builder.umi:
            write_tsv(os.path.join(t, "site_counts_prededup.tsv"), ["sample_id", "level"] + COUNT_COLS,
                      [[fmt(x) for x in r] for r in prededup_rows])
        rows = io.StringIO()
        w = csv.writer(rows, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "read_id", "class", "orientation", "molecules", "umis", "umi_1mismatch", "fates"])
        w.writerows(read_rows)
        gz_write(os.path.join(t, "reads.tsv.gz"), rows.getvalue())
        cat_rows = []
        for tg in targets:
            for site in tg["sites"]:
                exp_cats = design["site_types"][site["type"]].get("expected", {})
                for a in exp["analyses"]:
                    name = analysis_name(a)
                    cat_rows.append([name, tg["name"], site["pos"], site["type"], exp_cats.get(a, "-"),
                                     int(a in exp_cats)])
        write_tsv(os.path.join(t, "categories.tsv"),
                  ["analysis", "chr", "pos", "site_type", "expected", "clear_cut"], cat_rows)

        with open(os.path.join(edir, "samplesheet.csv"), "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["sample_id", "fastq", "treat", "rep", "celltype", "vector"])
            for s in samples:
                w.writerow([s["sample_id"], s["fastq"], s["treat"], s["rep"], s["celltype"], s["vector"]])
        write_analyses(os.path.join(edir, "analyses.yaml"), exp)
        write_params(os.path.join(edir, "params.yaml"), exp, params_root)

    write_readme(outdir, design)
    write_manifest(outdir)


def write_readme(outdir, design):
    exps = design["experiments"]
    lines = [
        "# Synthetic test data (generated; do not edit)",
        "",
        "Generated by `tests/simulate/simulate.py` from `tests/simulate/design.yaml` (seed %d):" % design["seed"],
        "",
        "    python tests/simulate/simulate.py --design tests/simulate/design.yaml --outdir tests/data/synthetic",
        "",
        "Same design and seed give byte-identical files (`MANIFEST.md5`). See `tests/simulate/README.md`",
        "for the design, read classes and truth tables.",
        "",
        "| Experiment | library_type | orientation_adapters | ont_adapter_mode | umi | samples |",
        "|---|---|---|---|---|---|",
    ]
    for name in sorted(exps):
        p, sm = exps[name]["params"], exps[name]["samples"]
        n = len(sm["vector"]) * len(sm["treat"]) * len(sm["rep"])
        lines.append("| `%s` | %s | %s | %s | %s | %d |" % (name, p["library_type"], p["orientation_adapters"],
                                                          p["ont_adapter_mode"], str(p["umi"]).lower(), n))
    lines += ["", "Run one experiment from the repo root:",
              "`nextflow run . -profile docker -params-file tests/data/synthetic/<experiment>/params.yaml`",
              "(paths in `params.yaml` are relative to the repo root; FASTQ paths in the sample sheet are",
              "relative to the sample sheet, per the §5.1 contract)."]
    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def analysis_name(a):
    return {"treatment": "WT_mod", "factor": "PUS7_dep_WT_v_KD"}[a]


def write_analyses(path, exp):
    lines = ["# Synthetic analyses (plan §5.4); expected calls in truth/categories.tsv", "analyses:"]
    for a in exp["analyses"]:
        if a == "treatment":
            lines += ["  - name: WT_mod", "    type: treatment", "    subset: {vector: [WT]}",
                      "    random_effects: \"\"", "    sesoi: 0.05", "    fdr: 0.05"]
        else:
            lines += ["  - name: PUS7_dep_WT_v_KD", "    type: factor", "    factor: vector",
                      "    levels: [KD, WT]", "    random_effects: \"\"", "    direction: positive",
                      "    sesoi: 0.05", "    fdr: 0.05"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_params(path, exp, root):
    """Params file; paths are relative to the repo root, where nextflow is launched."""
    p = exp["params"]
    ref = exp["reference"]
    d = "%s/%s" % (root, exp["name"])
    lines = [
        "# From the repo root: nextflow run . -profile docker -params-file %s/params.yaml" % d,
        "input: %s/samplesheet.csv" % d,
        "fasta: %s/reference/%s.fa" % (root, ref),
        "bed: %s/reference/%s.sites.bed" % (root, ref),
        "bed_coordinates: bed0",
        "library_type: %s" % p["library_type"],
        "orientation_adapters: %s" % p["orientation_adapters"],
        "ont_adapter_mode: %s" % p["ont_adapter_mode"],
        "umi: %s" % ("true" if p["umi"] else "false"),
        "analyses: %s/analyses.yaml" % d,
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_manifest(outdir):
    entries = []
    for root, _, files in os.walk(outdir):
        for f in files:
            if f == "MANIFEST.md5":
                continue
            p = os.path.join(root, f)
            with open(p, "rb") as fh:
                entries.append("%s  %s" % (hashlib.md5(fh.read()).hexdigest(), os.path.relpath(p, outdir)))
    with open(os.path.join(outdir, "MANIFEST.md5"), "w") as fh:
        fh.write("\n".join(sorted(entries, key=lambda e: e.split("  ", 1)[1])) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--design", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--params-root", default="tests/data/synthetic",
                    help="repo-relative location of the data, used for paths written into params.yaml")
    a = ap.parse_args(argv)
    with open(a.design) as fh:
        design = yaml.safe_load(fh)
    simulate(design, a.outdir, a.params_root)


if __name__ == "__main__":
    main()
