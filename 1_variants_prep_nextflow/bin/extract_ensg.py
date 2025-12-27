#!/usr/bin/env python3
"""
extract ensg ids from a gene list file.

args:
  genelist (str): path to input list
  output (str): path to write unique ensg ids

returns:
  None
"""

import argparse
import gzip
import re
from pathlib import Path
from typing import Iterable, Iterator

ENSG_RE = re.compile(r"(ENSG\d+(?:\.\d+)?)")


def _iter_lines(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            yield line


def extract_ensg(lines: Iterable[str]) -> list[str]:
    seen = set()
    ordered: list[str] = []
    for line in lines:
        for match in ENSG_RE.findall(line):
            base = match.split(".", 1)[0]
            if base in seen:
                continue
            seen.add(base)
            ordered.append(base)
    return ordered


def write_list(items: Iterable[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for item in items:
            handle.write(f"{item}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="extract ENSG IDs from a gene list")
    parser.add_argument("--genelist", required=True, help="path to input gene list")
    parser.add_argument("--output", required=True, help="path to write ensg list")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genelist = Path(args.genelist)
    if not genelist.exists():
        raise FileNotFoundError(f"genelist not found: {genelist}")
    items = extract_ensg(_iter_lines(genelist))
    if not items:
        raise ValueError("no ENSG IDs found in genelist")
    write_list(items, Path(args.output))


if __name__ == "__main__":
    main()


