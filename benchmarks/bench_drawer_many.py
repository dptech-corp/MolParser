from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path

import rdkit
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molparser.utils import drawer


CAPTIONS = [
    "CCOC(=O)c1ccccc1",
    "*c1ccccc1<sep><a>0:R[1]</a>",
    "*CCOCC*<sep><d>0:<dum></d><d>6:<dum></d>|Sg:n|",
    "COCCNCCOC<sep><g>[4:3]:[5:6]:|Sg:n|</g>",
    "CC*CC<sep><a>2:CH2?3</a>",
    "CC*<sep><a>2:<id>[blue]</a>",
    "C=CCC(C(C)*)*<sep><a>6:R[2]</a><a>7:R[1]</a><v>0::[0:2]</v>",
]


def unique_captions(count: int) -> list[str]:
    """Build structure-unique captions without relying on a render cache."""

    alphabet = "CNO"
    canonical_to_source: dict[str, str] = {}
    width = 8
    value = 0
    while len(canonical_to_source) < count:
        if value >= len(alphabet) ** width:
            width += 1
            value = 0
        current = value
        atoms = []
        for _ in range(width):
            atoms.append(alphabet[current % len(alphabet)])
            current //= len(alphabet)
        value += 1
        source = "".join(atoms)
        mol = Chem.MolFromSmiles(source)
        if mol is None:
            continue
        canonical_to_source.setdefault(
            Chem.MolToSmiles(mol, canonical=True),
            source,
        )

    result = []
    for index, source in enumerate(canonical_to_source.values()):
        atom_count = len(source)
        mode = index % 4
        if mode == 0:
            result.append(
                f"*{source}*<sep><d>0:<dum></d>"
                f"<d>{atom_count + 1}:<dum></d>|Sg:{2 + index % 199}|"
            )
        elif mode == 1:
            result.append(
                f"{source}<sep><g>[2:1]:"
                f"[{atom_count - 3}:{atom_count - 2}]:"
                f"|Sg:{2 + index % 199}|</g>"
            )
        elif mode == 2:
            result.append(
                f"{source}<sep><v>0::[0:{atom_count // 2}]</v>"
            )
        else:
            result.append(
                f"{source}*<sep><a>{atom_count}:<id>[blue]</a>"
            )
    return result


def measure(captions: list[str], workers: int, rounds: int) -> dict[str, float]:
    config = drawer.DrawingConfig()
    drawer.draw_many(captions[: min(len(captions), 14)], config=config, workers=workers)
    samples = []
    for _ in range(rounds):
        started = time.perf_counter()
        drawer.draw_many(captions, config=config, workers=workers)
        samples.append(time.perf_counter() - started)
    median = statistics.median(samples)
    sorted_samples = sorted(samples)
    p95 = sorted_samples[max(0, int(len(sorted_samples) * 0.95) - 1)]
    return {
        "median_seconds": round(median, 6),
        "p95_seconds": round(p95, 6),
        "median_images_per_second": round(len(captions) / median, 2),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=1400)
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument(
        "--dataset",
        choices=("repeated", "unique"),
        default="unique",
    )
    args = parser.parse_args()
    if args.count < 1 or args.rounds < 1:
        parser.error("--count and --rounds must be positive")

    if args.dataset == "unique":
        captions = unique_captions(args.count)
    else:
        captions = [CAPTIONS[index % len(CAPTIONS)] for index in range(args.count)]
    result = {
        "python": platform.python_version(),
        "rdkit": rdkit.__version__,
        "count": len(captions),
        "rounds": args.rounds,
        "dataset": args.dataset,
        "unique_captions": len(set(captions)),
        "workers": {
            str(workers): measure(captions, workers, args.rounds)
            for workers in args.workers
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
