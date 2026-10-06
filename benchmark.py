"""
Simulation benchmark: how well does the pipeline recover KNOWN, randomly placed
variants? Everything here is simulated from the reference, so it measures the
software (alignment + normalisation + naming), not biological accuracy.

Usage:  python benchmark.py --reference data/NM_000518.5.gb
"""
import argparse
import random
from collections import defaultdict
from pathlib import Path

import variant_caller as vc


def simulate(seq, rng, n_var, spacing=25):
    sites = sorted(rng.sample(range(5, len(seq) - spacing, spacing), n_var))
    truth = []
    for pos in sites:
        kind = rng.choice(["SNP", "SNP", "insertion", "deletion"])
        if kind == "SNP":
            truth.append(vc.Variant(pos, seq[pos], rng.choice([b for b in "ACGT" if b != seq[pos]])))
        elif kind == "insertion":
            truth.append(vc.Variant(pos, "", "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))))
        else:
            truth.append(vc.Variant(pos, seq[pos:pos + rng.randint(1, 6)], ""))
    sample = seq
    for v in sorted(truth, key=lambda v: -v.pos):
        sample = vc.apply_variant(sample, v)
    return sample, {vc.shift_left(seq, v) for v in truth}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--replicates", type=int, default=500)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()

    ref = vc.load_reference(Path(a.reference))
    rng = random.Random(a.seed)
    tp, fn, fp = defaultdict(int), defaultdict(int), defaultdict(int)
    roundtrip = 0
    for _ in range(a.replicates):
        sample, truth = simulate(ref.seq, rng, rng.randint(1, 5))
        called = set(vc.call_variants(ref.seq, sample).variants)
        for v in truth:
            (tp if v in called else fn)[v.vtype] += 1
        for v in called - truth:
            fp[v.vtype] += 1
        rebuilt = ref.seq
        for v in sorted(called, key=lambda v: -v.pos):
            rebuilt = vc.apply_variant(rebuilt, v)
        roundtrip += rebuilt == sample

    lines = [f"Simulation benchmark on {ref.accession}: {a.replicates} replicates, seed {a.seed}",
             f"{'type':<10}{'truth':>7}{'found':>7}{'missed':>8}{'false+':>8}{'recall':>9}{'precision':>11}"]
    for t in ("SNP", "insertion", "deletion"):
        n = tp[t] + fn[t]
        rec = tp[t] / n if n else float("nan")
        prec = tp[t] / (tp[t] + fp[t]) if (tp[t] + fp[t]) else float("nan")
        lines.append(f"{t:<10}{n:>7}{tp[t]:>7}{fn[t]:>8}{fp[t]:>8}{rec:>9.3f}{prec:>11.3f}")
    lines.append(f"Calls rebuild the sample exactly in {roundtrip}/{a.replicates} replicates")
    lines.append("CAVEAT: simulated, single clean sequences with well-separated variants. "
                 "Real data are harder.")
    text = "\n".join(lines)
    print(text)
    Path("results").mkdir(exist_ok=True)
    Path("results/benchmark.txt").write_text(text + "\n")


if __name__ == "__main__":
    main()
