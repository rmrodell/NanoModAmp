"""pytest checks on the WP1 simulator and the committed synthetic data (plan §7 WP1)."""
import csv
import gzip
import json
import os
import re
from collections import defaultdict

import pytest
import yaml

import simulate

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DESIGN = os.path.join(HERE, "design.yaml")
DATA = os.path.join(REPO, "tests", "data", "synthetic")


def load_design():
    with open(DESIGN) as fh:
        return yaml.safe_load(fh)


def tsv(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def experiments():
    return sorted(load_design()["experiments"])


@pytest.fixture(scope="module")
def regenerated(tmp_path_factory):
    out = tmp_path_factory.mktemp("sim") / "synthetic"
    simulate.simulate(load_design(), str(out))
    return out


def test_deterministic(regenerated, tmp_path):
    again = tmp_path / "again"
    simulate.simulate(load_design(), str(again))
    assert (regenerated / "MANIFEST.md5").read_bytes() == (again / "MANIFEST.md5").read_bytes()


def test_committed_data_is_current(regenerated):
    """tests/data/synthetic must be exactly what the recorded command produces."""
    with open(os.path.join(DATA, "MANIFEST.md5"), "rb") as fh:
        assert fh.read() == (regenerated / "MANIFEST.md5").read_bytes(), \
            "regenerate: python tests/simulate/simulate.py --design tests/simulate/design.yaml --outdir tests/data/synthetic"


def test_size_limit():
    total = sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(DATA) for f in fs)
    assert total <= 5 * 1024 * 1024


def test_gzip_headers_have_no_timestamp():
    for root, _, files in os.walk(DATA):
        for f in files:
            if f.endswith(".gz"):
                with open(os.path.join(root, f), "rb") as fh:
                    head = fh.read(10)
                assert head[4:8] == b"\x00\x00\x00\x00", f


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line
    return seqs


@pytest.mark.parametrize("ref", ["endo", "mpra"])
def test_sites_have_unambiguous_indels(ref):
    seqs = read_fasta(os.path.join(DATA, "reference", ref + ".fa"))
    for line in open(os.path.join(DATA, "reference", ref + ".sites.bed")):
        chrom, start, end = line.split("\t")[:3]
        start, end = int(start), int(end)
        assert end == start + 1, "single-site BED rows are standard 0-based (bed0)"
        s = seqs[chrom]
        assert s[start] == "T"
        assert s[start - 1] != "T" and s[start + 1] != "T"


def test_near_duplicate_pair():
    seqs = read_fasta(os.path.join(DATA, "reference", "mpra.fa"))
    assert seqs["NEARDUP_B"][:-1] == seqs["NEARDUP_A"]


@pytest.mark.parametrize("exp", experiments())
def test_samplesheet_contract(exp):
    rows = list(csv.DictReader(open(os.path.join(DATA, exp, "samplesheet.csv"))))
    assert list(rows[0])[:4] == ["sample_id", "fastq", "treat", "rep"]
    ids = [r["sample_id"] for r in rows]
    assert len(ids) == len(set(ids))
    for r in rows:
        assert re.fullmatch(r"[A-Za-z0-9._-]+", r["sample_id"])
        assert r["treat"] in ("input", "BS")
        assert os.path.isfile(os.path.join(DATA, exp, r["fastq"]))


@pytest.mark.parametrize("exp", experiments())
def test_params_match_schema(exp):
    schema = json.load(open(os.path.join(REPO, "nextflow_schema.json")))
    props = {}
    for d in schema.get("$defs", schema.get("definitions", {})).values():
        props.update(d.get("properties", {}))
    params = yaml.safe_load(open(os.path.join(DATA, exp, "params.yaml")))
    for k, v in params.items():
        assert k in props, k
        if "enum" in props[k]:
            assert v in props[k]["enum"], (k, v)
    for k in ("input", "fasta", "bed", "analyses"):
        assert os.path.exists(os.path.join(REPO, params[k])), k


@pytest.mark.parametrize("exp", experiments())
def test_analyses_match_schema(exp):
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.load(open(os.path.join(REPO, "assets", "schema_analyses.json")))
    jsonschema.validate(yaml.safe_load(open(os.path.join(DATA, exp, "analyses.yaml"))), schema)


@pytest.mark.parametrize("exp", experiments())
def test_funnel_consistent_with_reads(exp):
    params = load_design()["experiments"][exp]["params"]
    funnel = defaultdict(dict)
    for r in tsv(os.path.join(DATA, exp, "truth", "funnel.tsv")):
        funnel[r["sample_id"]][r["step"]] = int(r["records"])
    reads = tsv(os.path.join(DATA, exp, "truth", "reads.tsv.gz"))
    by_sample = defaultdict(list)
    for r in reads:
        by_sample[r["sample_id"]].append(r)
    for s, f in funnel.items():
        assert f["raw"] == len(by_sample[s])
        fates = [x for r in by_sample[s] for x in r["fates"].split(";") if x != "dropped_trim"]
        assert f["primary"] == len(fates)
        assert f["mapq_filtered"] == fates.count("mapped")
        if "merged" in f:
            assert f["merged"] == f["trim1_sense"] + f["trim1_antisense"]
            assert f["antisense_rc"] == f["trim1_antisense"]
        steps = list(f)
        if params["umi"]:
            mols = {m for r in by_sample[s] if "mapped" in r["fates"].split(";")
                    for m, fate in zip(r["molecules"].split(";"), r["fates"].split(";")) if fate == "mapped"}
            assert f["dedup"] == len(mols)
        else:
            assert "dedup" not in steps and "umi_extracted" not in steps
        if params["ont_adapter_mode"] == "three_prime_only":
            assert "trim1_3prime" in steps and "trim1_sense" not in steps
        if params["orientation_adapters"] == "pool" or params["library_type"] == "endogenous":
            assert "trim2_pool" not in steps


@pytest.mark.parametrize("exp", experiments())
def test_umi_design(exp):
    design = load_design()
    reads = tsv(os.path.join(DATA, exp, "truth", "reads.tsv.gz"))
    if not design["experiments"][exp]["params"]["umi"]:
        return
    mol_umi = {}
    variants = []
    for r in reads:
        if r["class"] not in ("valid", "multimapper"):
            continue
        m, u = r["molecules"], r["umis"]
        if r["umi_1mismatch"] == "1":
            variants.append((r["sample_id"], m, u))
        else:
            mol_umi.setdefault((r["sample_id"], m), u)
    by_target = defaultdict(list)
    for (s, m), u in mol_umi.items():
        by_target[(s, m.split(":")[0])].append(u)
    ham = lambda a, b: sum(x != y for x, y in zip(a, b))
    for umis in by_target.values():
        for i in range(len(umis)):
            for j in range(i + 1, len(umis)):
                assert ham(umis[i], umis[j]) >= design["min_umi_hamming"]
    assert variants, "design should include 1-mismatch UMI copies"
    for s, m, u in variants:
        assert ham(u, mol_umi[(s, m)]) == 1


@pytest.mark.parametrize("exp", experiments())
def test_site_counts(exp):
    rows = tsv(os.path.join(DATA, exp, "truth", "site_counts.tsv"))
    cov = load_design()["pipeline"]["min_coverage"]
    for r in rows:
        parts = [int(r[c]) for c in ("A.count", "C.count", "G.count", "T.count", "Deletion.count", "Insertion.count")]
        assert int(r["totalReads"]) == sum(parts)                     # D9
        assert abs(float(r["delrate"]) - int(r["Deletion.count"]) / int(r["totalReads"])) < 1e-12
        assert r["passes_min_coverage"] == str(int(int(r["totalReads"]) > cov))
        assert r["ref"] == "T" and r["strand"] == "+" and r["kmer"][2] == "T"
    totals = {int(r["totalReads"]) for r in rows}
    if exp in ("endogenous", "mpra"):
        assert {20, 21} <= totals, "coverage boundary sites (20 excluded, 21 kept)"
        assert any(int(r["Insertion.count"]) > 0 for r in rows)
        assert not any(r["chr"].startswith("NEARDUP") for r in rows), "multimapper reads are MAPQ 0"


@pytest.mark.parametrize("exp", experiments())
def test_categories(exp):
    rows = tsv(os.path.join(DATA, exp, "truth", "categories.tsv"))
    assert any(r["clear_cut"] == "1" for r in rows)
    for r in rows:
        if r["clear_cut"] == "1":
            assert r["expected"] in ("Modified", "Unmodified", "significant", "not_significant")


def test_three_prime_reads_lack_5prime_adapter():
    s5 = load_design()["adapters"]["ont_sense_5p"]
    for f in os.listdir(os.path.join(DATA, "legacy_3prime", "fastq")):
        with gzip.open(os.path.join(DATA, "legacy_3prime", "fastq", f), "rt") as fh:
            assert s5 not in fh.read()


def test_pool_library_has_no_ont_adapters():
    a = load_design()["adapters"]
    for f in os.listdir(os.path.join(DATA, "legacy_pool_noumi", "fastq")):
        with gzip.open(os.path.join(DATA, "legacy_pool_noumi", "fastq", f), "rt") as fh:
            txt = fh.read()
        assert a["ont_sense_5p"] not in txt and a["ont_sense_3p"] not in txt
