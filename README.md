# HBB variant-calling mini-pipeline

A small, transparent Python pipeline that compares a **sample sequence** to a **reference sequence** and reports **SNPs, insertions and deletions** in the human **beta-globin gene (HBB)** — the gene behind sickle cell disease and beta-thalassemia.

> \*\*Status: educational portfolio project. Not validated for clinical or diagnostic use.\*\*

## What it does

1. Downloads the reference **RefSeq NM\_000518.5** (HBB mRNA, 147-aa protein) from NCBI (or reads a local file).
2. Aligns each sample sequence to the reference (global alignment, affine gap penalties, Biopython `PairwiseAligner`). Reverse-complemented samples are detected automatically.
3. Extracts SNPs, insertions and deletions from the alignment and **normalises** indels (left-aligned for VCF, 3'-shifted for HGVS naming).
4. Names each variant in HGVS `c.` notation (e.g. `c.20A>T`) and predicts the protein consequence (missense, nonsense, frameshift, synonymous, UTR).
5. Writes, per sample: a **TSV** table, a **VCF** file and a **plot**; plus one combined table.

## Data (all public)

|What|Source|
|-|-|
|Reference sequence|NCBI RefSeq **NM\_000518.5**, *Homo sapiens* hemoglobin subunit beta (HBB), mRNA|
|Disease variants used for demo samples|ClinVar / dbSNP / ITHANET records, listed below|

|Demo variant|HGVS (NM\_000518.5)|rsID|Predicted protein change (checked against ClinVar/literature)|
|-|-|-|-|
|HbS (sickle)|c.20A>T|rs334|p.Glu7Val (also written Glu6Val in legacy numbering)|
|HbC|c.19G>A|rs33930165|p.Glu7Lys|
|HbE|c.79G>A|rs33950507|p.Glu27Lys|
|Codon 39|c.118C>T|rs11549407|p.Gln40Ter|
|Codon 8/9 (+G)|c.27dup|rs35699606|p.Ser10Valfs\*14|
|Codon 41/42 (-CTTT)|c.126\_129del|rs80356821|p.Phe42Leufs\*19|
|Codon 3 (synonymous, benign)|c.9T>C|rs713040|p.His3=|

**Important:** the demo samples are **simulated** — I applied these documented variants to the reference sequence myself. They are not sequences from patients.

## Quick start (Windows, Command Prompt)

```bat
py -3.14 -m venv .venv
.venv\\Scripts\\activate
pip install -r requirements.txt
python -m pytest tests -q
python variant\_caller.py --demo --email YOUR\_EMAIL@example.com
python benchmark.py --reference data\\NM\_000518.5.gb
```

Run your own sequence(s): `python variant\_caller.py --sample my\_samples.fasta`
(add `--free-end-gaps` if your sample is partial, e.g. coding sequence only).

No internet / NCBI blocked? Download the FASTA manually from the NCBI page for NM\_000518.5
(*Send to → Complete Record → File → FASTA*), save it as `data\\NM\_000518.5.fasta`, and add
`--reference data\\NM\_000518.5.fasta`. The coding sequence is then inferred as the longest ORF.

## Outputs (`results/`)

`<sample>.variants.tsv` · `<sample>.variants.vcf` · `<sample>.variants.png` · `all\_samples.variants.tsv` · `benchmark.txt`

## Validation

* **Unit tests (21)** cover SNPs, indels, strand detection, ambiguous bases, UTR variants, terminal variants,
normalisation, VCF correctness and a randomised simulation. Test CDS translation is checked against the published 147-aa beta-globin protein.
* **Demo check:** each simulated sample must be recovered with exactly the expected HGVS name(s).
* **Simulation benchmark** (`benchmark.py`): random SNPs/insertions/deletions (1–6 nt) are added to the reference and re-called.
Results from my run:

```
Simulation benchmark on NM\_000518.5: 500 replicates, seed 2026

type        truth  found  missed  false+   recall  precision

SNP           708    708       0       0    1.000      1.000

insertion     383    383       0       0    1.000      1.000

deletion      367    367       0       0    1.000      1.000

Calls rebuild the sample exactly in 500/500 replicates

CAVEAT: simulated, single clean sequences with well-separated variants. Real data are harder.


```

In my development runs, well-separated variants were recovered at \~100% recall/precision, and recall dropped to roughly 98–99%
when variants were packed only 3–10 nt apart. **This is simulated data on a short, clean sequence, so these numbers
overstate real-world performance.**

## Honest limitations

* **No real patient/sample data were analysed.** The demo is simulation-based, so a perfect score shows the code works, not that it is accurate on real biological data.
* **Sequence-vs-sequence, not read-based.** There are no sequencing reads, base qualities, depth, allele fractions or genotypes. Heterozygous (e.g. sickle-cell *trait*) samples cannot be represented by one sequence. Real pipelines (BWA + bcftools/GATK) do this properly.
* **The reference is an mRNA.** Introns, promoter and splice-site variants (e.g. IVS-I-1 `c.92+1G>A`, IVS-I-110 `c.93-21G>A`) and large gene deletions are invisible, and splicing effects are not predicted.
* **Closely spaced variants** (within \~10 nt) can occasionally be merged, split or named differently; normalisation is done against the reference only.
* **Terminal indels** are hidden when `--free-end-gaps` is used; without it, a partial sample shows its missing ends as large "deletions".
* **Protein prediction is basic:** SNPs and frameshifts only get full `p.` names. In-frame indels, start-loss and stop-loss get a consequence label without a `p.` name. Nonsense-mediated decay, splicing and protein stability are not modelled.
* **HGVS names are not validated with an official tool.** Cross-check important names with Mutalyzer or VariantValidator.
* **The "known variant" table is tiny and hand-curated.** It is not a clinical database and gives no pathogenicity classification; always consult ClinVar.
* **Alignment scores** (match +2, mismatch −3, gap open −6, extend −2) were chosen by hand, not tuned.
* **IUPAC ambiguity codes** (R, Y, …) are treated as `N` and ignored.
* **Test coverage gaps:** unit tests use the real HBB coding sequence with *synthetic* UTRs (random placeholder sequence). The live NCBI download function was not exercised in my development environment (NCBI was not reachable); I confirmed it only via the manual-FASTA route logic.

## Ideas for next steps

* Use a **genomic** HBB reference (RefSeqGene NG\_059281.1) so intron variants become visible, and compare against real GenBank isolate sequences.
* Add a read-based version (public SRA reads → BWA → bcftools) on WSL/Linux and compare to this sequence-based approach.
* Validate names with Mutalyzer; add zygosity-aware input.

## Files

`variant\_caller.py` (pipeline) · `benchmark.py` (simulation) · `tests/` (pytest) · `requirements.txt` · `.gitignore`

*Tools: Python, Biopython, NumPy, Matplotlib, pytest.*

