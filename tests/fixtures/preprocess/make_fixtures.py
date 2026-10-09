#!/usr/bin/env python3
"""Tiny hand-designed preprocessing fixtures for WP2 (plan §7 WP2 tests 1-14).

Deterministic: same bytes on every run (seeded sequences, gzip mtime 0). Run from this
directory:  python3 make_fixtures.py
Writes ref.fa, filter_test.sam, fastq/*.fastq.gz, fastq/S_dir/part{1,2}.fastq.gz,
empty_dir/.keep and expected.json (the truth the tests assert against).

Read layouts (sense orientation; antisense reads are the reverse complement of the whole read):
  ONT, MPRA       S5 + P5 + oligo + P3 + UMI + S3
  ONT, endogenous S5 + transcript + UMI + S3
  pool (legacy)   P5 + oligo + P3                  (no ONT adapters, no UMI; D29)
  3'-only legacy  transcript + UMI + S3            (no 5' ONT adapter; D31)
"""
import gzip
import json
import os
import random

S5, S3 = "TTTCTGTTGGTGCTGATATTGCG", "GAAGATAGAGCGACAGGCAAGT"
P5, P3 = "GACGCTCTTCCGATCT", "CACTCGGGCACCAAGGAC"
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


rng = random.Random(20261009)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def umi(existing):
    # 10-mer at Hamming distance >= 4 from every other UMI used
    while True:
        u = seq(10)
        if all(sum(a != b for a, b in zip(u, e)) >= 4 for e in existing):
            existing.append(u)
            return u


O1, O2 = seq(150), seq(150)
O2DUP = O2  # identical copy -> reads on O2 are multi-mappers (MAPQ 0)
E1 = seq(320)
ref = {"O1": O1, "O2": O2, "O2dup": O2DUP, "E1": E1}

umis = []
u1, u2, u3, u4, u5, uE1, uE2 = (umi(umis) for _ in range(7))
u1m = u1[:4] + ("A" if u1[4] != "A" else "C") + u1[5:]  # 1-mismatch copy of u1


def ont_mpra(ins, u):
    return S5 + P5 + ins + P3 + u + S3


def ont_endo(ins, u):
    return S5 + ins + u + S3


E1_INS = E1[40:260]  # 220 nt insert inside the transcript

samples = {
    # library_type mpra, defaults (ont, linked, umi)
    "mpra_ont": [
        ("dup1", ont_mpra(O1, u1)), ("dup2", ont_mpra(O1, u1)), ("dup3", ont_mpra(O1, u1)),
        ("dup_1mm", ont_mpra(O1, u1m)),
        ("sense_u2", ont_mpra(O1, u2)),
        ("antisense_u3", rc(ont_mpra(O1, u3))),
        ("only_5p", S5 + P5 + O1 + P3 + u4),
        ("only_3p", P5 + O1 + P3 + u4 + S3),
        ("too_short", ont_mpra(O1[:40], u4)),
        ("multimap_O2", ont_mpra(O2, u4)),
        ("chimera", ont_mpra(O1, u4) + rc(ont_mpra(O1, u5))),
    ],
    # library_type endogenous, defaults
    "endo_ont": [
        ("e_dup1", ont_endo(E1_INS, uE1)), ("e_dup2", ont_endo(E1_INS, uE1)),
        ("e_u2", ont_endo(E1_INS, uE2)),
        ("e_antisense_u2", rc(ont_endo(E1_INS, uE2))),
    ],
    # library_type mpra, orientation_adapters=pool, umi=false (D29, opt-in)
    "mpra_pool": [
        ("p_sense1", P5 + O1 + P3), ("p_sense2", P5 + O1 + P3),
        ("p_antisense", rc(P5 + O1 + P3)),
        ("p_one_adapter", P5 + O1),
    ],
    # library_type endogenous, ont_adapter_mode=three_prime_only (D31, opt-in)
    "endo_3prime": [
        ("t_3p_u1", E1_INS + uE1 + S3), ("t_3p_u2", E1_INS + uE2 + S3),
        ("t_no_adapter", E1[0:260]),
        ("t_too_short", E1_INS[:80] + uE1 + S3),
    ],
}


def fq(records):
    return "".join("@%s\n%s\n+\n%s\n" % (n, s, "I" * len(s)) for n, s in records)


def write_gz(path, text):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "wb") as fh:
        with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0, filename="") as gz:
            gz.write(text.encode())


with open("ref.fa", "w") as fh:
    for k, v in ref.items():
        fh.write(">%s\n%s\n" % (k, v))
for name, recs in samples.items():
    write_gz("fastq/%s.fastq.gz" % name, fq(recs))
# directory input: lexical order part1, part2 (test 8)
write_gz("fastq/S_dir/part2.fastq.gz", fq([("dir_b", ont_mpra(O1, u2))]))
write_gz("fastq/S_dir/part1.fastq.gz", fq([("dir_a", ont_mpra(O1, u1))]))
os.makedirs("empty_dir", exist_ok=True)
open("empty_dir/.keep", "w").close()

# SAM for the filter test (test 7): flags/MAPQ chosen around -F 2304 -q 30
q = "I" * 20
sam = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:O1\tLN:150"]
for name, flag, mapq in [("keep_mapq60", 0, 60), ("keep_mapq30", 0, 30), ("drop_mapq29", 0, 29),
                         ("drop_secondary", 256, 60), ("drop_supplementary", 2048, 60)]:
    sam.append("\t".join([name, str(flag), "O1", "11", str(mapq), "20M", "*", "0", "0", O1[10:30], q]))
with open("filter_test.sam", "w") as fh:
    fh.write("\n".join(sam) + "\n")

expected = {
    "inserts": {"O1": O1, "E1_INS": E1_INS},
    "umis": {"u1": u1, "u2": u2, "u3": u3, "uE1": uE1, "uE2": uE2},
    # read-funnel records per step (§5.2); steps a mode skips are absent
    "funnel": {
        "mpra_ont": {"raw": 11, "trim1_sense": 7, "trim1_antisense": 2, "antisense_rc": 2, "merged": 9,
                     "umi_extracted": 9, "trim2_pool": 9, "primary": 9, "mapq_filtered": 8, "dedup": 5},
        "endo_ont": {"raw": 4, "trim1_sense": 3, "trim1_antisense": 1, "antisense_rc": 1, "merged": 4,
                     "umi_extracted": 4, "primary": 4, "mapq_filtered": 4, "dedup": 2},
        "mpra_pool": {"raw": 4, "trim1_sense": 2, "trim1_antisense": 1, "antisense_rc": 1, "merged": 3,
                      "primary": 3, "mapq_filtered": 3},
        "endo_3prime": {"raw": 4, "trim1_3prime": 3, "umi_extracted": 3, "primary": 3, "mapq_filtered": 3,
                        "dedup": 3},
    },
    "both_orientations": {"mpra_ont": 1, "endo_ont": 0, "mpra_pool": 0},
    "filter_test_kept": ["keep_mapq30", "keep_mapq60"],
    "dir_order": ["dir_a", "dir_b"],
}
with open("expected.json", "w") as fh:
    json.dump(expected, fh, indent=1, sort_keys=True)
    fh.write("\n")
