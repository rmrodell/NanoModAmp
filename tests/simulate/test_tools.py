"""Check the simulator's mapping assumptions with minimap2 (skipped when minimap2 is absent).

For each experiment's first sample, every trimmed record that would reach the aligner is mapped
with the pipeline's default preset (D5). Expected:
  mapped   -> primary alignment to its target at position 1, MAPQ >= 30, no clipping, and exactly
              the designed indels (deletion of the site T, insertion right after it)
  mapq0    -> MAPQ < 30 (near-duplicate oligo pair)
  unmapped -> unmapped
Set MINIMAP2=/path/to/minimap2 if it is not on PATH.
"""
import os
import re
import shutil
import subprocess

import pytest
import yaml

import simulate

HERE = os.path.dirname(os.path.abspath(__file__))
MINIMAP2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
PRESET = {"endogenous": ["-ax", "splice", "-uf"], "mpra": ["-ax", "sr"]}

pytestmark = pytest.mark.skipif(not MINIMAP2, reason="minimap2 not available")


def design():
    with open(os.path.join(HERE, "design.yaml")) as fh:
        return yaml.safe_load(fh)


def indels(pos, cigar):
    """Reference positions of deletions and of bases followed by an insertion; clipping flag."""
    dels, ins, clipped = set(), set(), False
    ref = pos
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            ref += n
        elif op == "D":
            dels.update(range(ref, ref + n))
            ref += n
        elif op == "I":
            ins.add(ref - 1)
        elif op in "SH":
            clipped = True
        elif op == "N":
            dels.update(range(ref, ref + n))
            ref += n
    return dels, ins, clipped


@pytest.mark.parametrize("exp", sorted(design()["experiments"]))
def test_mapping_assumptions(exp, tmp_path):
    d = design()
    targets, params, recs = simulate.mapping_inputs(d, exp)
    ref = tmp_path / "ref.fa"
    ref.write_text("".join(">%s\n%s\n" % (t["name"], t["seq"]) for t in targets))
    fq = tmp_path / "reads.fa"
    fq.write_text("".join(">%s\n%s\n" % (rid, seq) for rid, seq, _, _ in recs))
    sam = subprocess.run([MINIMAP2] + PRESET[params["library_type"]] + [str(ref), str(fq)],
                         check=True, capture_output=True, text=True).stdout
    expect = {rid: (fate, mol) for rid, _, fate, mol in recs}
    seen = set()
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag, rname, pos, mapq, cigar = int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
        if flag & 0x900:
            continue
        fate, mol = expect[f[0]]
        seen.add(f[0])
        if fate == "unmapped":
            assert flag & 4, f[0]
            continue
        assert not flag & 4, f[0]
        if fate == "mapq0":
            assert mapq < 30, (f[0], mapq)
            continue
        assert rname == mol["target"] and pos == 1 and mapq >= 30, (f[0], rname, pos, mapq)
        dels, ins, clipped = indels(pos, cigar)
        assert not clipped, (f[0], cigar)
        assert dels == {p for p, o in mol["outcomes"].items() if o == "del"}, (f[0], cigar)
        assert ins == {p for p, o in mol["outcomes"].items() if o == "ins"}, (f[0], cigar)
    assert seen == set(expect)
