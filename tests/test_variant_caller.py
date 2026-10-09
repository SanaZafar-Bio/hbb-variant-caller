import random
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fixture import HBB_CDS, HBB_PROTEIN, make_fixture  # noqa: E402
import variant_caller as vc  # noqa: E402


@pytest.fixture(scope="module")
def ref(tmp_path_factory):
    p = make_fixture(tmp_path_factory.mktemp("ref") / "ref.gb")
    return vc.load_reference(p)


def apply_all(seq, variants):
    for v in sorted(variants, key=lambda v: -v.pos):
        seq = vc.apply_variant(seq, v)
    return seq


def names(ref, sample, **kw):
    res = vc.call_variants(ref.seq, sample, **kw)
    return [vc.hgvs_c(ref, v) for v in res.variants], res


# ---- reference / fixture sanity -------------------------------------------
def test_cds_translates_to_published_protein():
    assert str(Seq(HBB_CDS).translate()).rstrip("*") == HBB_PROTEIN


def test_reference_loads_cds_from_annotation(ref):
    assert (ref.cds_start, ref.cds_end, ref.cds_source) == (50, 494, "annotation")
    assert ref.c_pos(50) == "1" and ref.c_pos(49) == "-1" and ref.c_pos(494) == "*1"


def test_fasta_reference_uses_longest_orf(tmp_path, ref):
    f = tmp_path / "r.fasta"
    f.write_text(f">x\n{ref.seq}\n")
    r2 = vc.load_reference(f)
    assert r2.cds_source.startswith("longest ORF")
    assert r2.cds_start == 50 and r2.cds_end == 494


# ---- calling ---------------------------------------------------------------
def test_identical_sample_has_no_variants(ref):
    n, res = names(ref, ref.seq)
    assert n == [] and res.identity == 1.0


def test_snp_hbs(ref):
    s = vc.build_demo_sample(ref, [("sub", 20, "A", "T")])
    n, _ = names(ref, s)
    assert n == ["c.20A>T"]
    row = vc.annotate(ref, "t", vc.call_variants(ref.seq, s).variants)[0]
    assert row["hgvs_p"] == "p.Glu7Val" and row["known_rsid"] == "rs334"


def test_nonsense(ref):
    s = vc.build_demo_sample(ref, [("sub", 118, "C", "T")])
    row = vc.annotate(ref, "t", vc.call_variants(ref.seq, s).variants)[0]
    assert (row["consequence"], row["hgvs_p"]) == ("stop_gained", "p.Gln40Ter")


def test_deletion_matches_literature_name_and_protein(ref):
    s = vc.build_demo_sample(ref, [("del", 126, "CTTT", "")])
    res = vc.call_variants(ref.seq, s)
    row = vc.annotate(ref, "t", res.variants)[0]
    assert row["hgvs_c"] == "c.126_129del"
    assert row["hgvs_p"] == "p.Phe42Leufs*19"          # ClinVar: p.Phe42LeufsTer19


def test_insertion_is_named_as_duplication(ref):
    s = vc.build_demo_sample(ref, [("ins", 27, "", "G")])
    res = vc.call_variants(ref.seq, s)
    row = vc.annotate(ref, "t", res.variants)[0]
    assert row["hgvs_c"] == "c.27dup" and row["hgvs_p"] == "p.Ser10Valfs*14"


def test_non_duplication_insertion(ref):
    i = ref.c_to_index(60)  # insert 'CC' after c.60; check it's not a dup first
    ins = "CC" if ref.seq[i - 1:i + 1] != "CC" and ref.seq[i + 1:i + 3] != "CC" else "AA"
    s = ref.seq[:i + 1] + ins + ref.seq[i + 1:]
    n, _ = names(ref, s)
    assert len(n) == 1 and "ins" in n[0] or "dup" in n[0]
    assert apply_all(ref.seq, vc.call_variants(ref.seq, s).variants) == s


def test_utr_variants_are_labelled(ref):
    s = ref.seq[:10] + ("A" if ref.seq[10] != "A" else "C") + ref.seq[11:]
    row = vc.annotate(ref, "t", vc.call_variants(ref.seq, s).variants)[0]
    assert row["consequence"] == "5_prime_UTR_variant" and row["hgvs_c"].startswith("c.-")
    s = ref.seq[:-5] + ("A" if ref.seq[-5] != "A" else "C") + ref.seq[-4:]
    row = vc.annotate(ref, "t", vc.call_variants(ref.seq, s).variants)[0]
    assert row["consequence"] == "3_prime_UTR_variant" and row["hgvs_c"].startswith("c.*")


def test_reverse_complement_detected(ref):
    s = vc.reverse_complement(vc.build_demo_sample(ref, [("sub", 20, "A", "T")]))
    n, res = names(ref, s)
    assert n == ["c.20A>T"] and res.reverse_complemented


def test_ambiguous_bases_are_not_called(ref):
    s = ref.seq[:200] + "N" + ref.seq[201:]
    res = vc.call_variants(ref.seq, s)
    assert res.variants == [] and res.n_ambiguous == 1


def test_variants_at_first_and_last_base(ref):
    flip = lambda b: "A" if b != "A" else "C"
    s = flip(ref.seq[0]) + ref.seq[1:-1] + flip(ref.seq[-1])
    res = vc.call_variants(ref.seq, s)
    assert [v.pos for v in res.variants] == [0, len(ref.seq) - 1]


def test_multi_variant_roundtrip(ref):
    s = vc.build_demo_sample(ref, [("sub", 9, "T", "C"), ("sub", 79, "G", "A"),
                                   ("del", 126, "CTTT", "")])
    res = vc.call_variants(ref.seq, s)
    assert len(res.variants) == 3
    assert apply_all(ref.seq, res.variants) == s


def test_partial_sample_needs_free_end_gaps(ref):
    cds_only = ref.seq[ref.cds_start:ref.cds_end]
    strict = vc.call_variants(ref.seq, cds_only)
    assert any(v.vtype == "deletion" for v in strict.variants)       # flanks look like deletions
    free = vc.call_variants(ref.seq, cds_only, free_end_gaps=True)
    assert free.variants == [] and "flank" in free.flank_note


# ---- normalisation ---------------------------------------------------------
def test_shift_left_right_homopolymer():
    seq = "GGTCAAAAGCTT"
    d = vc.Variant(5, "A", "")                    # delete one A in AAAA (indices 4-7)
    left, right = vc.shift_left(seq, d), vc.shift_right(seq, d)
    assert left.pos == 4 and right.pos == 7
    assert vc.apply_variant(seq, left) == vc.apply_variant(seq, right)
    ins = vc.Variant(6, "", "A")
    l2, r2 = vc.shift_left(seq, ins), vc.shift_right(seq, ins)
    assert l2.pos == 4 and r2.pos == 8
    assert vc.apply_variant(seq, l2) == vc.apply_variant(seq, r2)


# ---- output files ----------------------------------------------------------
def test_vcf_roundtrip_reproduces_sample(ref, tmp_path):
    s = vc.build_demo_sample(ref, [("sub", 9, "T", "C"), ("ins", 27, "", "G"),
                                   ("del", 126, "CTTT", "")])
    res = vc.call_variants(ref.seq, s)
    out = tmp_path / "x.vcf"
    vc.write_vcf(out, ref, "t", res.variants)
    recs = [l.split("\t") for l in out.read_text().splitlines() if not l.startswith("#")]
    assert len(recs) == 3
    seq = ref.seq
    for chrom, pos, _id, r, a, *_ in sorted(recs, key=lambda x: -int(x[1])):
        p = int(pos) - 1
        assert seq[p:p + len(r)] == r, "VCF REF does not match the reference"
        seq = seq[:p] + a + seq[p + len(r):]
    assert seq == s


def test_tsv_written(ref, tmp_path):
    rows = vc.annotate(ref, "t", vc.call_variants(
        ref.seq, vc.build_demo_sample(ref, [("sub", 20, "A", "T")])).variants)
    out = tmp_path / "x.tsv"
    vc.write_tsv(out, rows)
    assert out.read_text().splitlines()[0].split("\t") == vc.TSV_FIELDS


def test_demo_builder_rejects_wrong_reference(ref):
    with pytest.raises(ValueError):
        vc.build_demo_sample(ref, [("sub", 20, "G", "T")])  # c.20 is A, not G


def test_end_to_end_demo(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    refpath = make_fixture(tmp_path / "ref.gb")
    code = vc.main(["--reference", str(refpath), "--demo", "--outdir", str(tmp_path / "out")])
    assert code == 0
    assert (tmp_path / "out" / "all_samples.variants.tsv").exists()
    assert any(p.suffix == ".png" for p in (tmp_path / "out").iterdir())


# ---- randomised simulation -------------------------------------------------
def simulate(ref, rng, n_var):
    seq = ref.seq
    L = len(seq)
    sites = sorted(rng.sample(range(5, L - 25, 25), n_var))
    truth = []
    for pos in sites:
        kind = rng.choice(["snp", "snp", "ins", "del"])
        if kind == "snp":
            alt = rng.choice([b for b in "ACGT" if b != seq[pos]])
            truth.append(vc.Variant(pos, seq[pos], alt))
        elif kind == "ins":
            truth.append(vc.Variant(pos, "", "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))))
        else:
            truth.append(vc.Variant(pos, seq[pos:pos + rng.randint(1, 6)], ""))
    sample = apply_all(seq, truth)
    return sample, {vc.shift_left(seq, v) for v in truth}


def test_random_simulation_recovers_variants(ref):
    rng = random.Random(2026)
    total = found = exact_roundtrip = 0
    reps = 150
    for _ in range(reps):
        sample, truth = simulate(ref, rng, rng.randint(1, 5))
        res = vc.call_variants(ref.seq, sample)
        called = set(res.variants)
        total += len(truth)
        found += len(truth & called)
        exact_roundtrip += apply_all(ref.seq, res.variants) == sample
    assert exact_roundtrip == reps               # calls always rebuild the sample exactly
    assert found / total >= 0.97                 # >=97% exact match to simulated truth
