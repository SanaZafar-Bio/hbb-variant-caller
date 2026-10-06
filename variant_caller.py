"""
HBB variant-calling mini-pipeline
=================================
Compares a SAMPLE sequence to a REFERENCE sequence (default: human HBB mRNA,
RefSeq NM_000518.5), and reports SNPs, insertions and deletions.

Method (deliberately simple and transparent):
  1. Global pairwise alignment (Biopython PairwiseAligner, affine gaps).
  2. Walk the alignment: mismatches -> SNPs, gaps -> insertions / deletions.
  3. Normalise indels (left-aligned for VCF, 3'-shifted for HGVS, as each standard requires).
  4. Name variants in HGVS "c." notation and predict the protein consequence.
  5. Write TSV + VCF + a plot.

This is an educational pipeline for single, clean sequences. It is NOT a
replacement for read-based callers (GATK, bcftools, DeepVariant) and must
not be used for clinical decisions. See README.md for limitations.
"""
from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqUtils import seq3

REF_ACCESSION = "NM_000518.5"

# Tiny hand-curated lookup (NOT a database). Keys are HGVS c. names on NM_000518.5.
# rsIDs/names were checked against ClinVar / ITHANET / dbSNP pages in Oct 2026.
KNOWN_VARIANTS = {
    "c.20A>T": ("rs334", "HbS (sickle cell variant)"),
    "c.19G>A": ("rs33930165", "HbC"),
    "c.79G>A": ("rs33950507", "HbE"),
    "c.118C>T": ("rs11549407", "codon 39 (C>T) nonsense, beta-thalassemia"),
    "c.27dup": ("rs35699606", "codon 8/9 (+G) frameshift, beta-thalassemia"),
    "c.126_129del": ("rs80356821", "codon 41/42 (-CTTT) frameshift, beta-thalassemia"),
    "c.9T>C": ("rs713040", "codon 3 (T>C) synonymous, benign"),
}

# ----------------------------------------------------------------------------
# Reference handling
# ----------------------------------------------------------------------------
@dataclass
class Reference:
    accession: str
    seq: str
    cds_start: int  # 0-based, first base of ATG
    cds_end: int    # 0-based exclusive, includes the stop codon
    cds_source: str  # "annotation" or "longest ORF (inferred)"

    def c_pos(self, pos0: int) -> str:
        """0-based transcript index -> HGVS c. position string."""
        if pos0 < self.cds_start:
            return f"-{self.cds_start - pos0}"
        if pos0 >= self.cds_end:
            return f"*{pos0 - self.cds_end + 1}"
        return str(pos0 - self.cds_start + 1)

    def c_to_index(self, c: int) -> int:
        """Positive c. coordinate (inside CDS) -> 0-based transcript index."""
        return self.cds_start + c - 1


def longest_orf(seq: str) -> tuple[int, int]:
    """Longest forward ATG...stop ORF. Returns (start, end_exclusive)."""
    stops = {"TAA", "TAG", "TGA"}
    best = (0, 0)
    for start in range(len(seq) - 2):
        if seq[start:start + 3] != "ATG":
            continue
        for i in range(start, len(seq) - 2, 3):
            if seq[i:i + 3] in stops:
                if i + 3 - start > best[1] - best[0]:
                    best = (start, i + 3)
                break
    if best == (0, 0):
        raise ValueError("No ATG...stop open reading frame found in reference.")
    return best


def load_reference(path: Path) -> Reference:
    """Load a GenBank (.gb) or FASTA (.fasta) reference. GenBank CDS is preferred."""
    suffix = path.suffix.lower()
    fmt = "genbank" if suffix in (".gb", ".gbk", ".genbank") else "fasta"
    rec = SeqIO.read(str(path), fmt)
    seq = clean_sequence(str(rec.seq))
    if fmt == "genbank":
        cds = [f for f in rec.features if f.type == "CDS"]
        if cds:
            loc = cds[0].location
            return Reference(rec.id, seq, int(loc.start), int(loc.end), "annotation")
    start, end = longest_orf(seq)
    return Reference(rec.id, seq, start, end, "longest ORF (inferred)")


def fetch_reference(accession: str, email: str, outdir: Path) -> Path:
    """Download a GenBank record from NCBI (cached on disk)."""
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{accession}.gb"
    if out.exists() and out.stat().st_size > 0:
        print(f"[ref] using cached {out}")
        return out
    Entrez.email = email
    print(f"[ref] downloading {accession} from NCBI ...")
    with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as h:
        text = h.read()
    if not text.startswith("LOCUS"):
        raise RuntimeError(f"NCBI did not return a GenBank record for {accession}:\n{text[:200]}")
    out.write_text(text)
    return out


def clean_sequence(s: str) -> str:
    """Uppercase, U->T, anything not A/C/G/T -> N."""
    s = s.upper().replace("U", "T").replace(" ", "").replace("\n", "")
    return "".join(c if c in "ACGT" else "N" for c in s)


# ----------------------------------------------------------------------------
# Variant representation + normalisation
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class Variant:
    pos: int   # 0-based. For insertions: index of the reference base the insert sits BEFORE.
    ref: str   # "" for pure insertion
    alt: str   # "" for pure deletion

    @property
    def vtype(self) -> str:
        if self.ref and self.alt:
            if len(self.ref) == len(self.alt) == 1:
                return "SNP"
            return "MNP" if len(self.ref) == len(self.alt) else "complex"
        return "insertion" if self.alt else "deletion"


def shift_left(seq: str, v: Variant) -> Variant:
    """Left-align a pure indel (VCF convention)."""
    pos, ref, alt = v.pos, v.ref, v.alt
    if not ref and alt:
        while pos > 0 and seq[pos - 1] == alt[-1]:
            alt = alt[-1] + alt[:-1]
            pos -= 1
    elif ref and not alt:
        while pos > 0 and seq[pos - 1] == ref[-1]:
            ref = ref[-1] + ref[:-1]
            pos -= 1
    return Variant(pos, ref, alt)


def shift_right(seq: str, v: Variant) -> Variant:
    """3'-most placement of a pure indel (HGVS convention)."""
    pos, ref, alt = v.pos, v.ref, v.alt
    if not ref and alt:
        while pos < len(seq) and seq[pos] == alt[0]:
            alt = alt[1:] + alt[0]
            pos += 1
    elif ref and not alt:
        end = pos + len(ref)
        while end < len(seq) and seq[pos] == seq[end]:
            ref = ref[1:] + seq[end]
            pos += 1
            end += 1
    return Variant(pos, ref, alt)


# ----------------------------------------------------------------------------
# Alignment + variant calling
# ----------------------------------------------------------------------------
def make_aligner(free_end_gaps: bool = False) -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2
    a.mismatch_score = -3
    a.open_gap_score = -6
    a.extend_gap_score = -2
    if free_end_gaps:  # for partial sequences (e.g. CDS-only sample vs full mRNA)
        a.end_gap_score = 0
    return a


def reverse_complement(s: str) -> str:
    return str(Seq(s).reverse_complement())


@dataclass
class CallResult:
    variants: list
    n_ambiguous: int
    identity: float
    aligned_ref_bases: int
    flank_note: str
    reverse_complemented: bool


def call_variants(ref: str, sample: str, free_end_gaps: bool = False) -> CallResult:
    aligner = make_aligner(free_end_gaps)
    fwd = aligner.score(ref, sample)
    rev_sample = reverse_complement(sample)
    rc = aligner.score(ref, rev_sample) > fwd
    if rc:
        sample = rev_sample

    aln = aligner.align(ref, sample)[0]
    tblocks, qblocks = aln.aligned
    tb = [(int(a), int(b)) for a, b in tblocks]
    qb = [(int(a), int(b)) for a, b in qblocks]

    raw: list[Variant] = []
    n_amb = 0
    matches = 0
    aligned_cols = 0
    for (rs, re_), (qs, qe) in zip(tb, qb):
        for i in range(re_ - rs):
            r, q = ref[rs + i], sample[qs + i]
            aligned_cols += 1
            if r == q and r != "N":
                matches += 1
            elif q == "N" or r == "N":
                n_amb += 1
            else:
                raw.append(Variant(rs + i, r, q))

    # Gaps between aligned blocks (and, unless free_end_gaps, at the two ends)
    flank_note = ""
    t_edges = [(0, 0)] + tb + [(len(ref), len(ref))]
    q_edges = [(0, 0)] + qb + [(len(sample), len(sample))]
    for k in range(len(t_edges) - 1):
        r0, r1 = t_edges[k][1], t_edges[k + 1][0]
        q0, q1 = q_edges[k][1], q_edges[k + 1][0]
        rseq, qseq = ref[r0:r1], sample[q0:q1]
        if not rseq and not qseq:
            continue
        is_end = k == 0 or k == len(t_edges) - 2
        if free_end_gaps and is_end:
            flank_note += (f" unaligned {'5-prime' if k == 0 else '3-prime'} flank "
                           f"(ref {len(rseq)} nt, sample {len(qseq)} nt);")
            continue
        if "N" in qseq:
            n_amb += 1
            continue
        raw.append(Variant(r0, rseq, qseq))

    # Normalise (left-aligned canonical form), de-duplicate, sort
    norm = sorted({shift_left(ref, v) for v in raw}, key=lambda v: (v.pos, v.ref, v.alt))
    identity = matches / aligned_cols if aligned_cols else 0.0
    return CallResult(norm, n_amb, identity, aligned_cols, flank_note.strip(), rc)


# ----------------------------------------------------------------------------
# HGVS naming + protein consequence
# ----------------------------------------------------------------------------
def apply_variant(seq: str, v: Variant) -> str:
    return seq[:v.pos] + v.alt + seq[v.pos + len(v.ref):]


def hgvs_c(refobj: Reference, v: Variant) -> str:
    r = shift_right(refobj.seq, v)  # HGVS 3' rule
    cp = refobj.c_pos
    if r.vtype == "SNP":
        return f"c.{cp(r.pos)}{r.ref}>{r.alt}"
    if r.vtype == "deletion":
        a, b = cp(r.pos), cp(r.pos + len(r.ref) - 1)
        return f"c.{a}del" if a == b else f"c.{a}_{b}del"
    if r.vtype == "insertion":
        n = len(r.alt)
        if r.pos >= n and refobj.seq[r.pos - n:r.pos] == r.alt:  # duplication
            a, b = cp(r.pos - n), cp(r.pos - 1)
            return f"c.{a}dup" if a == b else f"c.{a}_{b}dup"
        return f"c.{cp(r.pos - 1)}_{cp(r.pos)}ins{r.alt}"
    a, b = cp(r.pos), cp(r.pos + len(r.ref) - 1)
    span = a if a == b else f"{a}_{b}"
    return f"c.{span}delins{r.alt}"


def _translate_to_stop(nt: str) -> str:
    nt = nt[: len(nt) // 3 * 3]
    aa = str(Seq(nt).translate())
    i = aa.find("*")
    return aa if i == -1 else aa[: i + 1]


def predict_protein(refobj: Reference, v: Variant) -> tuple[str, str]:
    """Return (consequence, HGVS p.) for a variant. p. is '' when not computed."""
    s, e = v.pos, v.pos + len(v.ref)  # affected reference interval
    if e <= refobj.cds_start and not (v.vtype == "insertion" and s == refobj.cds_start):
        return "5_prime_UTR_variant", ""
    if s >= refobj.cds_end:
        return "3_prime_UTR_variant", ""
    if s < refobj.cds_start:
        return "start_region_complex", ""

    ref_aa = _translate_to_stop(refobj.seq[refobj.cds_start:refobj.cds_end])
    mutated = apply_variant(refobj.seq, v)
    alt_aa = _translate_to_stop(mutated[refobj.cds_start:])
    if not alt_aa.startswith("M"):
        return "start_lost", ""

    delta = len(v.alt) - len(v.ref)
    first = next((i for i in range(min(len(ref_aa), len(alt_aa)))
                  if ref_aa[i] != alt_aa[i]), None)
    if first is None and len(ref_aa) != len(alt_aa):
        first = min(len(ref_aa), len(alt_aa))

    if delta % 3 != 0:  # frameshift
        i = first if first is not None else 0
        ra = seq3(ref_aa[i]) if i < len(ref_aa) else "Ter"
        aa_alt = alt_aa[i] if i < len(alt_aa) else "*"
        ext = f"*{len(alt_aa) - i}" if alt_aa.endswith("*") else "*?"
        return "frameshift_variant", f"p.{ra}{i + 1}{seq3(aa_alt)}fs{ext}"

    if v.vtype == "SNP":
        codon_no = (s - refobj.cds_start) // 3 + 1
        if first is None:
            return "synonymous_variant", f"p.{seq3(ref_aa[codon_no - 1])}{codon_no}="
        ra, qa = ref_aa[first], alt_aa[first] if first < len(alt_aa) else "*"
        if qa == "*":
            return "stop_gained", f"p.{seq3(ra)}{first + 1}Ter"
        if ra == "*":
            return "stop_lost", ""
        return "missense_variant", f"p.{seq3(ra)}{first + 1}{seq3(qa)}"

    kind = {"deletion": "inframe_deletion", "insertion": "inframe_insertion"}.get(
        v.vtype, "inframe_delins")
    return kind, ""


def context(refobj: Reference, v: Variant, flank: int = 8) -> str:
    left = refobj.seq[max(0, v.pos - flank):v.pos]
    right = refobj.seq[v.pos + len(v.ref):v.pos + len(v.ref) + flank]
    return f"{left}[{v.ref or '-'}>{v.alt or '-'}]{right}"


def annotate(refobj: Reference, sample_id: str, variants: list) -> list[dict]:
    rows = []
    for v in variants:
        name = hgvs_c(refobj, v)
        csq, p = predict_protein(refobj, v)
        rsid, label = KNOWN_VARIANTS.get(name, (".", ""))
        rows.append({
            "sample": sample_id,
            "type": v.vtype,
            "pos_1based_vcf": v.pos + 1,
            "ref": v.ref or "-",
            "alt": v.alt or "-",
            "hgvs_c": name,
            "consequence": csq,
            "hgvs_p": p,
            "known_rsid": rsid,
            "known_label": label,
            "context": context(refobj, v),
        })
    return rows


# ----------------------------------------------------------------------------
# Output writers
# ----------------------------------------------------------------------------
TSV_FIELDS = ["sample", "type", "pos_1based_vcf", "ref", "alt", "hgvs_c", "consequence",
              "hgvs_p", "known_rsid", "known_label", "context"]


def write_tsv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TSV_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _info_escape(s: str) -> str:
    return s.replace("=", "%3D").replace(";", "%3B").replace(",", "%2C").replace(" ", "_")


def write_vcf(path: Path, refobj: Reference, sample_id: str, variants: list) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        f"##source=variant_caller.py ({sample_id})",
        f"##contig=<ID={refobj.accession},length={len(refobj.seq)}>",
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">',
        '##INFO=<ID=HGVSC,Number=1,Type=String,Description="HGVS c. name (3-prime rule)">',
        '##INFO=<ID=CSQ,Number=1,Type=String,Description="Predicted consequence">',
        '##INFO=<ID=HGVSP,Number=1,Type=String,Description="Predicted protein change">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for row, v in zip(annotate(refobj, sample_id, variants), variants):
        seq = refobj.seq
        if v.vtype in ("SNP", "MNP", "complex"):
            pos, ref, alt = v.pos + 1, v.ref, v.alt
        elif v.pos > 0:  # indel with preceding anchor base
            anchor = seq[v.pos - 1]
            pos, ref, alt = v.pos, anchor + v.ref, anchor + v.alt
        else:            # indel at position 1: anchor is the following base
            nxt = seq[v.pos + len(v.ref)]
            pos, ref, alt = 1, v.ref + nxt, v.alt + nxt
        info = f"TYPE={row['type']};HGVSC={_info_escape(row['hgvs_c'])};CSQ={row['consequence']}"
        if row["hgvs_p"]:
            info += f";HGVSP={_info_escape(row['hgvs_p'])}"
        rsid = row["known_rsid"]
        lines.append(f"{refobj.accession}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t{info}")
    path.write_text("\n".join(lines) + "\n")


def plot_variants(path: Path, refobj: Reference, sample_id: str, rows: list[dict],
                  variants: list) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"SNP": "#d62728", "insertion": "#2ca02c", "deletion": "#1f77b4",
              "MNP": "#9467bd", "complex": "#ff7f0e"}
    L = len(refobj.seq)
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.add_patch(plt.Rectangle((0, -0.12), L, 0.24, fc="#e8e8e8", ec="none"))
    ax.add_patch(plt.Rectangle((refobj.cds_start, -0.2), refobj.cds_end - refobj.cds_start,
                               0.4, fc="#9ecae1", ec="#3182bd"))
    ax.text((refobj.cds_start + refobj.cds_end) / 2, 0, "CDS (147 aa)" if
            (refobj.cds_end - refobj.cds_start) == 444 else "CDS",
            ha="center", va="center", fontsize=9)
    for i, (row, v) in enumerate(zip(rows, variants)):
        x = v.pos + (0.5 if v.ref else 0)
        h = 0.9 + 0.5 * (i % 3)
        ax.vlines(x, 0.2, h, color=colors[row["type"]], lw=1.5)
        ax.plot(x, h, "o", color=colors[row["type"]], ms=7)
        ax.text(x, h + 0.08, row["hgvs_c"], rotation=35, ha="left", va="bottom", fontsize=8)
    for t, c in colors.items():
        ax.plot([], [], "o", color=c, label=t)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=8, ncol=5,
              frameon=False)
    ax.set_xlim(0, L)
    ax.set_ylim(-0.4, 2.6)
    ax.set_yticks([])
    ax.set_xlabel(f"Position on {refobj.accession} (nt)")
    ax.set_title(f"Variants called in sample '{sample_id}'  (n={len(variants)})")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Demo samples (simulated from documented ClinVar variants)
# ----------------------------------------------------------------------------
# Each edit: (kind, c_position, ref_bases_expected, alt_bases). c. = CDS coordinate.
DEMO_SAMPLES = {
    "HbS_c.20A>T":        [("sub", 20, "A", "T")],
    "HbC_c.19G>A":        [("sub", 19, "G", "A")],
    "HbE_c.79G>A":        [("sub", 79, "G", "A")],
    "CD39_c.118C>T":      [("sub", 118, "C", "T")],
    "CD8-9_c.27dup":      [("ins", 27, "", "G")],       # insert G after c.27
    "CD41-42_c.126_129del": [("del", 126, "CTTT", "")],
    "multi_c.9T>C+c.79G>A+c.126_129del": [("sub", 9, "T", "C"), ("sub", 79, "G", "A"),
                                          ("del", 126, "CTTT", "")],
}


def build_demo_sample(refobj: Reference, edits: list) -> str:
    seq = refobj.seq
    for kind, c, exp, alt in sorted(edits, key=lambda e: -e[1]):  # right-to-left
        i = refobj.c_to_index(c)
        if kind in ("sub", "del"):
            found = seq[i:i + len(exp)]
            if found != exp:
                raise ValueError(
                    f"Reference mismatch at c.{c}: expected {exp}, found {found}. "
                    "Is the reference really NM_000518.5?")
            seq = seq[:i] + alt + seq[i + len(exp):]
        else:  # insertion after c. position c
            seq = seq[:i + 1] + alt + seq[i + 1:]
    return seq


def write_demo_fasta(refobj: Reference, path: Path) -> dict:
    truth = {}
    with open(path, "w") as fh:
        for name, edits in DEMO_SAMPLES.items():
            s = build_demo_sample(refobj, edits)
            fh.write(f">{name}\n{s}\n")
            truth[name] = sorted(_expected_names(name))
        # reverse-complemented copy to exercise strand detection
        s = build_demo_sample(refobj, DEMO_SAMPLES["HbS_c.20A>T"])
        fh.write(f">HbS_reverse_complemented\n{reverse_complement(s)}\n")
        truth["HbS_reverse_complemented"] = ["c.20A>T"]
    return truth


def _expected_names(sample_name: str) -> list:
    """Expected HGVS names for a demo sample = the names embedded in the sample id."""
    core = sample_name.split("_", 1)[1]
    return core.split("+")


# ----------------------------------------------------------------------------
# Command-line interface
# ----------------------------------------------------------------------------
def run(args) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.reference:
        ref_path = Path(args.reference)
    else:
        if not args.email:
            sys.exit("ERROR: NCBI requires an email. Add:  --email you@example.com  "
                     "(or download the FASTA manually and pass --reference, see README).")
        ref_path = fetch_reference(REF_ACCESSION, args.email, Path("data"))
    refobj = load_reference(ref_path)
    print(f"[ref] {refobj.accession}: {len(refobj.seq)} nt, CDS {refobj.cds_start + 1}-"
          f"{refobj.cds_end} ({refobj.cds_source})")

    truth = None
    if args.demo:
        Path("data").mkdir(exist_ok=True)
        sample_path = Path("data") / "demo_samples.fasta"
        truth = write_demo_fasta(refobj, sample_path)
        print(f"[demo] wrote simulated samples -> {sample_path}")
    elif args.sample:
        sample_path = Path(args.sample)
    else:
        sys.exit("ERROR: give --sample FILE.fasta or use --demo")

    all_rows = []
    n_correct = 0
    for rec in SeqIO.parse(str(sample_path), "fasta"):
        sid = rec.id
        sample_seq = clean_sequence(str(rec.seq))
        res = call_variants(refobj.seq, sample_seq, args.free_end_gaps)
        rows = annotate(refobj, sid, res.variants)
        all_rows.extend(rows)
        safe = "".join(c if c.isalnum() or c in "-._" else "_" for c in sid)
        write_tsv(outdir / f"{safe}.variants.tsv", rows)
        write_vcf(outdir / f"{safe}.variants.vcf", refobj, sid, res.variants)
        plot_variants(outdir / f"{safe}.variants.png", refobj, sid, rows, res.variants)

        print(f"\n== {sid}: {len(sample_seq)} nt | identity over aligned columns "
              f"{res.identity:.3%} | {len(rows)} variant(s)"
              f"{' | reverse-complemented' if res.reverse_complemented else ''}")
        if res.n_ambiguous:
            print(f"   note: {res.n_ambiguous} ambiguous (N) position(s) ignored")
        if res.flank_note:
            print(f"   note: {res.flank_note}")
        for r in rows:
            extra = f"  [{r['known_rsid']} {r['known_label']}]" if r["known_rsid"] != "." else ""
            print(f"   {r['hgvs_c']:<16} {r['type']:<9} {r['consequence']:<20} "
                  f"{r['hgvs_p']:<18}{extra}")
        if truth is not None:
            got = sorted(r["hgvs_c"] for r in rows)
            ok = got == truth[sid]
            n_correct += ok
            print(f"   demo check: {'PASS' if ok else 'FAIL'}  expected {truth[sid]}")

    write_tsv(outdir / "all_samples.variants.tsv", all_rows)
    if truth is not None:
        print(f"\n[demo] {n_correct}/{len(truth)} simulated samples recovered exactly.")
        print("[demo] NOTE: samples are simulated from the reference, so this checks that the "
              "code works, not that it is accurate on real biological data.")
    print(f"\nResults written to {outdir.resolve()}")
    return 0 if truth is None or n_correct == len(truth) else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="HBB variant-calling mini-pipeline")
    p.add_argument("--reference", help="Reference .gb/.gbk or .fasta (default: download NM_000518.5)")
    p.add_argument("--email", help="Your email (required by NCBI when downloading)")
    p.add_argument("--sample", help="Sample FASTA (one or more records)")
    p.add_argument("--demo", action="store_true", help="Run on simulated ClinVar-based samples")
    p.add_argument("--outdir", default="results")
    p.add_argument("--free-end-gaps", action="store_true",
                   help="Do not call terminal gaps (use when sample is partial, e.g. CDS only)")
    return run(p.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
