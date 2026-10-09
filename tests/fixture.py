"""Test fixture: HBB CDS + synthetic UTRs written as a GenBank file.

IMPORTANT: the UTRs here are random placeholder sequence (seeded), NOT the real
NM_000518.5 UTRs. The CDS is the human beta-globin coding sequence; its
translation is checked against the published 147-aa protein in the tests.
"""
import random
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

HBB_CDS = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGG"
    "CTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAG"
    "GCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCAC"
    "TGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGCGTGCTGGCCCATCACTTTGGCAAAGAATTCACC"
    "CCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
)
HBB_PROTEIN = (  # 147 aa, as listed for NP_000509 / P68871
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAV"
    "MGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVL"
    "VCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
)


def make_fixture(path: Path, utr5=50, utr3=134, seed=7) -> Path:
    rng = random.Random(seed)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    seq = rnd(utr5) + HBB_CDS + rnd(utr3)
    rec = SeqRecord(Seq(seq), id="TEST_HBB.1", name="TEST_HBB",
                    description="fixture: HBB CDS + synthetic UTRs",
                    annotations={"molecule_type": "mRNA"})
    rec.features.append(SeqFeature(FeatureLocation(utr5, utr5 + len(HBB_CDS), strand=1), type="CDS"))
    SeqIO.write(rec, str(path), "genbank")
    return path
