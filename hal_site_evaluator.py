#!/usr/bin/env python3
"""Scalable, run-length HAL site-naming evaluation.

Input prediction/truth TSV columns:
    block_id  species  chrom  start  end  strand

Coordinates are zero-based, half-open.  The ``map`` command batches records by
species and calls halLiftover once per species.  The ``score`` command compares
the resulting ancestral-coordinate runs without enumerating individual bases.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

FIELDS = ("block_id", "species", "chrom", "start", "end", "strand")
RUN_FIELDS = (
    "block_id", "species", "ancestor", "anc_chrom", "anc_start", "anc_end", "strand"
)


def read_tsv(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    return rows


def write_tsv(path: Path, fields: tuple[str, ...], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_source(rows: list[dict[str, str]]) -> None:
    seen = set()
    for line_no, row in enumerate(rows, 2):
        key = (row["species"], row["block_id"])
        if key in seen:
            raise ValueError(
                f"line {line_no}: duplicate (species, block_id) {key}; "
                "merge the occurrence first or give it a unique block_id"
            )
        seen.add(key)
        start, end = int(row["start"]), int(row["end"])
        if start < 0 or end <= start:
            raise ValueError(f"line {line_no}: invalid interval [{start}, {end})")
        if row["strand"] not in {"+", "-", "."}:
            raise ValueError(f"line {line_no}: strand must be +, -, or .")


def map_to_ancestor(input_tsv: Path, hal: Path, ancestor: str, output_tsv: Path,
                    hal_liftover: str) -> None:
    rows = read_tsv(input_tsv, FIELDS)
    validate_source(rows)
    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_species[row["species"]].append(row)

    if shutil.which(hal_liftover) is None:
        raise RuntimeError(f"cannot find executable: {hal_liftover}")

    mapped = []
    with tempfile.TemporaryDirectory(prefix="hal-site-eval-") as tmp:
        tmpdir = Path(tmp)
        for species, species_rows in sorted(by_species.items()):
            source_bed = tmpdir / f"{species}.source.bed"
            target_bed = tmpdir / f"{species}.ancestor.bed"
            # BED name carries the block ID through split mappings.
            with source_bed.open("w") as handle:
                for row in species_rows:
                    handle.write(
                        f'{row["chrom"]}\t{row["start"]}\t{row["end"]}\t'
                        f'{row["block_id"]}\t0\t{row["strand"]}\n'
                    )
            cmd = [hal_liftover, str(hal), species, str(source_bed), ancestor, str(target_bed)]
            completed = subprocess.run(cmd, text=True, capture_output=True)
            if completed.returncode:
                raise RuntimeError(
                    f"halLiftover failed for {species} (exit {completed.returncode}):\n"
                    f"{completed.stderr.strip()}"
                )
            if not target_bed.exists():
                continue
            with target_bed.open() as handle:
                for line_no, line in enumerate(handle, 1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 4:
                        raise ValueError(f"{target_bed}:{line_no}: expected at least 4 BED columns")
                    mapped.append({
                        "block_id": cols[3], "species": species, "ancestor": ancestor,
                        "anc_chrom": cols[0], "anc_start": cols[1], "anc_end": cols[2],
                        "strand": cols[5] if len(cols) >= 6 else ".",
                    })
    write_tsv(output_tsv, RUN_FIELDS, mapped)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[list[int]] = []
    for start, end in sorted(intervals):
        if not result or start > result[-1][1]:
            result.append([start, end])
        else:
            result[-1][1] = max(result[-1][1], end)
    return [(a, b) for a, b in result]


def interval_length(by_chrom: dict[str, list[tuple[int, int]]]) -> int:
    return sum(end - start for xs in by_chrom.values() for start, end in merge_intervals(xs))


def intersection_length(a: dict[str, list[tuple[int, int]]],
                        b: dict[str, list[tuple[int, int]]]) -> int:
    total = 0
    for chrom in a.keys() & b.keys():
        xs, ys = merge_intervals(a[chrom]), merge_intervals(b[chrom])
        i = j = 0
        while i < len(xs) and j < len(ys):
            total += max(0, min(xs[i][1], ys[j][1]) - max(xs[i][0], ys[j][0]))
            if xs[i][1] <= ys[j][1]:
                i += 1
            else:
                j += 1
    return total


def load_runs(path: Path):
    rows = read_tsv(path, RUN_FIELDS)
    runs = defaultdict(lambda: defaultdict(list))
    strands = defaultdict(list)
    for row in rows:
        key = (row["species"], row["block_id"])
        start, end = int(row["anc_start"]), int(row["anc_end"])
        if end <= start:
            raise ValueError(f"{path}: invalid mapped interval [{start}, {end})")
        runs[key][row["anc_chrom"]].append((start, end))
        if row["strand"] in {"+", "-"}:
            strands[key].append(row["strand"])
    return runs, strands


def dominant_strand(values: list[str]) -> str:
    if not values:
        return "."
    plus, minus = values.count("+"), values.count("-")
    return "+" if plus > minus else "-" if minus > plus else "."


def score_runs(pred_path: Path, truth_path: Path, output_tsv: Path) -> dict[str, float | int]:
    pred, pred_strands = load_runs(pred_path)
    truth, truth_strands = load_runs(truth_path)
    keys = sorted(pred.keys() | truth.keys())
    output = []
    total_intersection = total_pred = total_truth = 0
    orientation_correct = orientation_total = 0
    for species, block_id in keys:
        key = (species, block_id)
        p, t = pred.get(key, {}), truth.get(key, {})
        p_len, t_len = interval_length(p), interval_length(t)
        inter = intersection_length(p, t)
        union = p_len + t_len - inter
        precision = inter / p_len if p_len else 0.0
        recall = inter / t_len if t_len else 0.0
        jaccard = inter / union if union else 1.0
        ps, ts = dominant_strand(pred_strands[key]), dominant_strand(truth_strands[key])
        orient = "NA"
        if ps != "." and ts != ".":
            orient = "1" if ps == ts else "0"
            orientation_correct += int(ps == ts)
            orientation_total += 1
        output.append({
            "block_id": block_id, "species": species, "pred_bp": p_len,
            "truth_bp": t_len, "intersection_bp": inter, "union_bp": union,
            "precision": f"{precision:.6f}", "recall": f"{recall:.6f}",
            "jaccard": f"{jaccard:.6f}", "orientation_correct": orient,
        })
        total_intersection += inter
        total_pred += p_len
        total_truth += t_len
    fields = ("block_id", "species", "pred_bp", "truth_bp", "intersection_bp", "union_bp",
              "precision", "recall", "jaccard", "orientation_correct")
    write_tsv(output_tsv, fields, output)
    total_union = total_pred + total_truth - total_intersection
    return {
        "n_occurrences": len(keys),
        "micro_precision": total_intersection / total_pred if total_pred else 0.0,
        "micro_recall": total_intersection / total_truth if total_truth else 0.0,
        "micro_jaccard": total_intersection / total_union if total_union else 1.0,
        "orientation_accuracy": (
            orientation_correct / orientation_total if orientation_total else None
        ),
        "orientation_n": orientation_total,
    }


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="hal-site-test-") as tmp:
        root = Path(tmp)
        pred = root / "pred.tsv"
        truth = root / "truth.tsv"
        out = root / "scores.tsv"
        # B1: P=[120,200)+[250,280)=110; T=[100,200)+[250,300)=150.
        pred_rows = [
            {"block_id": "B1", "species": "eagle", "ancestor": "avianAnc",
             "anc_chrom": "Anc4", "anc_start": "120", "anc_end": "200", "strand": "-"},
            {"block_id": "B1", "species": "eagle", "ancestor": "avianAnc",
             "anc_chrom": "Anc4", "anc_start": "250", "anc_end": "280", "strand": "-"},
        ]
        truth_rows = [
            {"block_id": "B1", "species": "eagle", "ancestor": "avianAnc",
             "anc_chrom": "Anc4", "anc_start": "100", "anc_end": "200", "strand": "-"},
            {"block_id": "B1", "species": "eagle", "ancestor": "avianAnc",
             "anc_chrom": "Anc4", "anc_start": "250", "anc_end": "300", "strand": "-"},
        ]
        write_tsv(pred, RUN_FIELDS, pred_rows)
        write_tsv(truth, RUN_FIELDS, truth_rows)
        summary = score_runs(pred, truth, out)
        assert abs(summary["micro_precision"] - 1.0) < 1e-12
        assert abs(summary["micro_recall"] - 110 / 150) < 1e-12
        assert abs(summary["micro_jaccard"] - 110 / 150) < 1e-12
        assert summary["orientation_accuracy"] == 1.0
    print("SELF-TEST PASSED")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    m = sub.add_parser("map", help="map block intervals to one HAL ancestral genome")
    m.add_argument("--input", type=Path, required=True)
    m.add_argument("--hal", type=Path, required=True)
    m.add_argument("--ancestor", required=True)
    m.add_argument("--output", type=Path, required=True)
    m.add_argument("--hal-liftover", default="halLiftover")
    s = sub.add_parser("score", help="compare two mapped run tables")
    s.add_argument("--pred", type=Path, required=True)
    s.add_argument("--truth", type=Path, required=True)
    s.add_argument("--output", type=Path, required=True)
    sub.add_parser("self-test", help="test interval scoring without HAL")
    return p


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "map":
            map_to_ancestor(args.input, args.hal, args.ancestor, args.output, args.hal_liftover)
        elif args.command == "score":
            print(json.dumps(score_runs(args.pred, args.truth, args.output), indent=2))
        else:
            self_test()
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
