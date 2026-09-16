#!/usr/bin/env python3
"""Run-length HAL ancestral-coordinate evaluation for synteny blocks.

Coordinates are zero-based and half-open throughout.  The ``map`` command
batches records by species and calls ``halLiftover`` once per species.  The
``score`` command preserves the legacy same-label comparison, while ``match``
compares independently named block systems by ancestral interval overlap.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable

__version__ = "0.2.0"

FIELDS = ("block_id", "species", "chrom", "start", "end", "strand")
SOURCE_FIELDS = (
    "block_id", "occurrence_id", "species", "chrom", "start", "end", "strand",
    "copy_id", "status",
)
RUN_FIELDS = (
    "block_id", "species", "ancestor", "anc_chrom", "anc_start", "anc_end", "strand"
)
MAPPED_FIELDS = (
    "block_id", "occurrence_id", "species", "ancestor", "anc_chrom", "anc_start",
    "anc_end", "strand", "copy_id", "mapping_status", "source_status",
    "source_chrom", "source_start", "source_end", "source_strand",
    "source_length", "mapped_bp", "unique_bp", "duplicated_bp",
    "ambiguous_bp", "unmapped_bp",
)
STATUSES = {"unique", "duplicated", "ambiguous", "unmapped"}


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    strand: str = "."
    species: str = ""
    occurrence_id: str = ""
    status: str = "unique"

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Occurrence:
    block_id: str
    occurrence_id: str
    species: str
    chrom: str
    start: int
    end: int
    strand: str
    copy_id: str = "1"
    status: str = "unspecified"

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Family:
    block_id: str
    intervals: list[Interval] = field(default_factory=list)
    occurrence_ids: set[str] = field(default_factory=set)
    species: set[str] = field(default_factory=set)
    occurrence_status: dict[str, str] = field(default_factory=dict)


def read_tsv(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def write_tsv(path: Path, fields: tuple[str, ...], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _stable_occurrence_id(row: dict[str, str], ordinal: int) -> str:
    payload = "\t".join(
        [row["block_id"], row["species"], row["chrom"], row["start"], row["end"], row["strand"], str(ordinal)]
    )
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]
    return f'{row["block_id"]}:{row["species"]}:{digest}'


def parse_int(value: str, label: str, line_no: int | None = None) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        prefix = f"line {line_no}: " if line_no else ""
        raise ValueError(f"{prefix}{label} must be an integer") from exc


def validate_source(rows: list[dict[str, str]]) -> list[Occurrence]:
    occurrences: list[Occurrence] = []
    seen_occurrences: set[str] = set()
    generated_counts: dict[tuple[str, str, str, str, str, str], int] = defaultdict(int)
    for line_no, row in enumerate(rows, 2):
        for field_name in FIELDS:
            if not row.get(field_name, ""):
                raise ValueError(f"line {line_no}: missing required value {field_name}")
        start = parse_int(row["start"], "start", line_no)
        end = parse_int(row["end"], "end", line_no)
        if start < 0 or end <= start:
            raise ValueError(f"line {line_no}: invalid interval [{start}, {end})")
        if row["strand"] not in {"+", "-", "."}:
            raise ValueError(f"line {line_no}: strand must be +, -, or .")
        base_key = (
            row["block_id"], row["species"], row["chrom"], row["start"], row["end"], row["strand"]
        )
        generated_counts[base_key] += 1
        occurrence_id = (row.get("occurrence_id") or "").strip()
        if not occurrence_id:
            occurrence_id = _stable_occurrence_id(row, generated_counts[base_key])
        if occurrence_id in seen_occurrences:
            raise ValueError(f"line {line_no}: duplicate occurrence_id {occurrence_id!r}")
        seen_occurrences.add(occurrence_id)
        occurrences.append(Occurrence(
            block_id=row["block_id"], occurrence_id=occurrence_id, species=row["species"],
            chrom=row["chrom"], start=start, end=end, strand=row["strand"],
            copy_id=(row.get("copy_id") or "1"), status=(row.get("status") or "unspecified"),
        ))
    return occurrences


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
        total += intersection_length_lists(a[chrom], b[chrom])
    return total


def intersection_length_lists(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    xs, ys = merge_intervals(a), merge_intervals(b)
    total = i = j = 0
    while i < len(xs) and j < len(ys):
        total += max(0, min(xs[i][1], ys[j][1]) - max(xs[i][0], ys[j][0]))
        if xs[i][1] <= ys[j][1]:
            i += 1
        else:
            j += 1
    return total


def intervals_by_chrom(intervals: Iterable[Interval]) -> dict[str, list[tuple[int, int]]]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for interval in intervals:
        by_chrom[interval.chrom].append((interval.start, interval.end))
    return by_chrom


def classify_liftover_status(source: Occurrence, runs: list[Interval]) -> tuple[str, dict[str, int]]:
    source_len = source.length
    if not runs:
        return "unmapped", {
            "mapped_bp": 0, "unique_bp": 0, "duplicated_bp": 0,
            "ambiguous_bp": 0, "unmapped_bp": source_len,
        }
    mapped_bp = sum(run.length for run in runs)
    chroms = {run.chrom for run in runs}
    strands = {run.strand for run in runs if run.strand in {"+", "-"}}
    merged_bp = interval_length(intervals_by_chrom(runs))
    if len(chroms) > 1 or len(strands) > 1:
        status = "ambiguous"
    elif mapped_bp > source_len or mapped_bp > merged_bp:
        status = "duplicated"
    elif mapped_bp < source_len:
        status = "ambiguous"
    else:
        status = "unique"
    stats = {
        "mapped_bp": mapped_bp,
        "unique_bp": mapped_bp if status == "unique" else 0,
        "duplicated_bp": mapped_bp if status == "duplicated" else 0,
        "ambiguous_bp": mapped_bp if status == "ambiguous" else 0,
        "unmapped_bp": max(0, source_len - min(source_len, mapped_bp)),
    }
    return status, stats


def map_to_ancestor(input_tsv: Path, hal: Path, ancestor: str, output_tsv: Path,
                    hal_liftover: str) -> None:
    rows = read_tsv(input_tsv, FIELDS)
    occurrences = validate_source(rows)
    by_species: dict[str, list[Occurrence]] = defaultdict(list)
    for occurrence in occurrences:
        by_species[occurrence.species].append(occurrence)

    if shutil.which(hal_liftover) is None:
        raise RuntimeError(f"cannot find executable: {hal_liftover}")

    mapped_rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="hal-site-eval-") as tmp:
        tmpdir = Path(tmp)
        for species, species_occurrences in sorted(by_species.items()):
            source_bed = tmpdir / f"{species}.source.bed"
            target_bed = tmpdir / f"{species}.ancestor.bed"
            record_to_occurrence = {
                f"HSE{idx:012d}": occurrence
                for idx, occurrence in enumerate(species_occurrences, 1)
            }
            with source_bed.open("w") as handle:
                for record_id, occurrence in record_to_occurrence.items():
                    handle.write(
                        f"{occurrence.chrom}\t{occurrence.start}\t{occurrence.end}\t"
                        f"{record_id}\t0\t{occurrence.strand}\n"
                    )
            cmd = [hal_liftover, str(hal), species, str(source_bed), ancestor, str(target_bed)]
            completed = subprocess.run(cmd, text=True, capture_output=True)
            if completed.returncode:
                raise RuntimeError(
                    f"halLiftover failed for {species} (exit {completed.returncode}):\n"
                    f"{completed.stderr.strip()}"
                )
            lifted: dict[str, list[Interval]] = defaultdict(list)
            if target_bed.exists():
                with target_bed.open() as handle:
                    for line_no, line in enumerate(handle, 1):
                        if not line.strip() or line.startswith("#"):
                            continue
                        cols = line.rstrip("\n").split("\t")
                        if len(cols) < 4:
                            raise ValueError(f"{target_bed}:{line_no}: expected at least 4 BED columns")
                        record_id = cols[3]
                        if record_id not in record_to_occurrence:
                            raise ValueError(f"{target_bed}:{line_no}: unknown lifted record {record_id!r}")
                        start, end = int(cols[1]), int(cols[2])
                        if end <= start:
                            raise ValueError(f"{target_bed}:{line_no}: invalid mapped interval [{start}, {end})")
                        occurrence = record_to_occurrence[record_id]
                        lifted[record_id].append(Interval(
                            chrom=cols[0], start=start, end=end,
                            strand=cols[5] if len(cols) >= 6 else ".",
                            species=species, occurrence_id=occurrence.occurrence_id,
                        ))
            for record_id, occurrence in record_to_occurrence.items():
                runs = sorted(lifted.get(record_id, []), key=lambda x: (x.chrom, x.start, x.end, x.strand))
                status, stats = classify_liftover_status(occurrence, runs)
                if not runs:
                    mapped_rows.append(_mapped_row(occurrence, ancestor, "", "", "", ".", status, stats))
                for run in runs:
                    mapped_rows.append(_mapped_row(
                        occurrence, ancestor, run.chrom, run.start, run.end, run.strand, status, stats
                    ))
    write_tsv(output_tsv, MAPPED_FIELDS, mapped_rows)


def _mapped_row(occurrence: Occurrence, ancestor: str, anc_chrom: str, anc_start: int | str,
                anc_end: int | str, strand: str, mapping_status: str,
                stats: dict[str, int]) -> dict[str, object]:
    return {
        "block_id": occurrence.block_id, "occurrence_id": occurrence.occurrence_id,
        "species": occurrence.species, "ancestor": ancestor, "anc_chrom": anc_chrom,
        "anc_start": anc_start, "anc_end": anc_end, "strand": strand,
        "copy_id": occurrence.copy_id, "mapping_status": mapping_status,
        "source_status": occurrence.status, "source_chrom": occurrence.chrom,
        "source_start": occurrence.start, "source_end": occurrence.end,
        "source_strand": occurrence.strand, "source_length": occurrence.length,
        **stats,
    }


def load_runs(path: Path):
    rows = read_tsv(path, RUN_FIELDS)
    runs = defaultdict(lambda: defaultdict(list))
    strands = defaultdict(list)
    for row in rows:
        if not row.get("anc_chrom"):
            continue
        key = (row["species"], row["block_id"])
        start, end = int(row["anc_start"]), int(row["anc_end"])
        if start < 0 or end <= start:
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


def score_runs(pred_path: Path, truth_path: Path, output_tsv: Path) -> dict[str, float | int | None]:
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
        "orientation_accuracy": orientation_correct / orientation_total if orientation_total else None,
        "orientation_n": orientation_total,
    }


def load_families(path: Path) -> tuple[dict[str, Family], str | None]:
    rows = read_tsv(path, RUN_FIELDS)
    families: dict[str, Family] = {}
    ancestor: str | None = None
    occurrence_seen: dict[str, str] = {}
    for line_no, row in enumerate(rows, 2):
        block_id = row["block_id"]
        species = row["species"]
        occurrence_id = row.get("occurrence_id") or _legacy_mapped_occurrence_id(row)
        status = row.get("mapping_status") or row.get("status") or "unique"
        if status not in STATUSES:
            status = "ambiguous"
        previous_block = occurrence_seen.setdefault(occurrence_id, block_id)
        if previous_block != block_id:
            raise ValueError(f"{path}: line {line_no}: occurrence_id {occurrence_id!r} spans multiple block_ids")
        family = families.setdefault(block_id, Family(block_id=block_id))
        family.occurrence_ids.add(occurrence_id)
        family.species.add(species)
        family.occurrence_status[occurrence_id] = status
        if row.get("ancestor"):
            ancestor = ancestor or row["ancestor"]
            if ancestor != row["ancestor"]:
                raise ValueError(f"{path}: multiple ancestors found ({ancestor!r}, {row['ancestor']!r})")
        if status == "unmapped" or not row.get("anc_chrom"):
            continue
        start = parse_int(row["anc_start"], "anc_start", line_no)
        end = parse_int(row["anc_end"], "anc_end", line_no)
        if start < 0 or end <= start:
            raise ValueError(f"{path}: line {line_no}: invalid mapped interval [{start}, {end})")
        family.intervals.append(Interval(
            chrom=row["anc_chrom"], start=start, end=end, strand=row.get("strand") or ".",
            species=species, occurrence_id=occurrence_id, status=status,
        ))
    return families, ancestor


def _legacy_mapped_occurrence_id(row: dict[str, str]) -> str:
    payload = "\t".join([
        row["block_id"], row["species"], row.get("ancestor", ""), row.get("anc_chrom", ""),
        row.get("anc_start", ""), row.get("anc_end", ""), row.get("strand", "."),
    ])
    return f'{row["block_id"]}:{row["species"]}:{hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]}'


def family_len(family: Family) -> int:
    return interval_length(intervals_by_chrom(family.intervals))


def overlap_details(left: Family, right: Family) -> dict[str, object]:
    inter = 0
    orientation_same = orientation_total = 0
    species: set[str] = set()
    supporting_occurrences: set[tuple[str, str]] = set()
    status_counts: Counter[str] = Counter()
    left_by_chrom: dict[str, list[Interval]] = defaultdict(list)
    right_by_chrom: dict[str, list[Interval]] = defaultdict(list)
    for interval in left.intervals:
        left_by_chrom[interval.chrom].append(interval)
    for interval in right.intervals:
        right_by_chrom[interval.chrom].append(interval)
    for chrom in sorted(left_by_chrom.keys() & right_by_chrom.keys()):
        for li in left_by_chrom[chrom]:
            for ri in right_by_chrom[chrom]:
                bp = max(0, min(li.end, ri.end) - max(li.start, ri.start))
                if bp == 0:
                    continue
                inter += bp
                if li.species == ri.species:
                    species.add(li.species)
                supporting_occurrences.add(("left", li.occurrence_id))
                supporting_occurrences.add(("right", ri.occurrence_id))
                status_counts[li.status] += 1
                status_counts[ri.status] += 1
                if li.strand in {"+", "-"} and ri.strand in {"+", "-"}:
                    orientation_total += bp
                    orientation_same += bp if li.strand == ri.strand else 0
    total_support = sum(status_counts.values())
    unique = status_counts["unique"]
    nonunique = status_counts["duplicated"] + status_counts["ambiguous"]
    return {
        "intersection_bp": inter,
        "supporting_species": len(species),
        "orientation_agreement": orientation_same / orientation_total if orientation_total else None,
        "fraction_unique_support": unique / total_support if total_support else 0.0,
        "fraction_duplicated_or_ambiguous": nonunique / total_support if total_support else 0.0,
    }


def candidate_edges(left: dict[str, Family], right: dict[str, Family],
                    thresholds: dict[str, float]) -> list[dict[str, object]]:
    left_len = {block_id: family_len(family) for block_id, family in left.items()}
    right_len = {block_id: family_len(family) for block_id, family in right.items()}
    by_chrom_left: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    by_chrom_right: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for block_id, family in left.items():
        for chrom, runs in intervals_by_chrom(family.intervals).items():
            for start, end in merge_intervals(runs):
                by_chrom_left[chrom].append((start, end, block_id))
    for block_id, family in right.items():
        for chrom, runs in intervals_by_chrom(family.intervals).items():
            for start, end in merge_intervals(runs):
                by_chrom_right[chrom].append((start, end, block_id))
    pair_seen: set[tuple[str, str]] = set()
    for chrom in sorted(by_chrom_left.keys() & by_chrom_right.keys()):
        left_runs = sorted(by_chrom_left[chrom])
        right_runs = sorted(by_chrom_right[chrom])
        active_start = 0
        for l_start, l_end, l_id in left_runs:
            while active_start < len(right_runs) and right_runs[active_start][1] <= l_start:
                active_start += 1
            j = active_start
            while j < len(right_runs) and right_runs[j][0] < l_end:
                r_start, r_end, r_id = right_runs[j]
                if min(l_end, r_end) > max(l_start, r_start):
                    pair_seen.add((l_id, r_id))
                j += 1
    edges: list[dict[str, object]] = []
    for left_id, right_id in sorted(pair_seen):
        details = overlap_details(left[left_id], right[right_id])
        inter = int(details["intersection_bp"])
        l_len, r_len = left_len[left_id], right_len[right_id]
        union = l_len + r_len - inter
        jaccard = inter / union if union else 1.0
        left_cov = inter / l_len if l_len else 0.0
        right_cov = inter / r_len if r_len else 0.0
        if inter < thresholds["min_intersection_bp"]:
            continue
        if jaccard < thresholds["min_jaccard"]:
            continue
        if left_cov < thresholds["min_left_coverage"]:
            continue
        if right_cov < thresholds["min_right_coverage"]:
            continue
        if int(details["supporting_species"]) < thresholds["min_supporting_species"]:
            continue
        edges.append({
            "left_block_id": left_id, "right_block_id": right_id,
            "left_bp": l_len, "right_bp": r_len, "intersection_bp": inter,
            "union_bp": union, "jaccard": f"{jaccard:.6f}",
            "left_coverage": f"{left_cov:.6f}", "right_coverage": f"{right_cov:.6f}",
            "orientation_agreement": _fmt_optional(details["orientation_agreement"]),
            "supporting_species": details["supporting_species"],
            "fraction_unique_support": f"{details['fraction_unique_support']:.6f}",
            "fraction_duplicated_or_ambiguous": f"{details['fraction_duplicated_or_ambiguous']:.6f}",
        })
    return edges


def _fmt_optional(value: object) -> str:
    return "NA" if value is None else f"{float(value):.6f}"


def connected_components(edges: list[dict[str, object]]) -> list[tuple[list[str], list[str], list[dict[str, object]]]]:
    left_to_right: dict[str, set[str]] = defaultdict(set)
    right_to_left: dict[str, set[str]] = defaultdict(set)
    edge_by_pair = {(e["left_block_id"], e["right_block_id"]): e for e in edges}
    for edge in edges:
        left_to_right[str(edge["left_block_id"])].add(str(edge["right_block_id"]))
        right_to_left[str(edge["right_block_id"])].add(str(edge["left_block_id"]))
    visited_left: set[str] = set()
    visited_right: set[str] = set()
    components = []
    for seed in sorted(left_to_right):
        if seed in visited_left:
            continue
        q = deque([("L", seed)])
        left_nodes: set[str] = set()
        right_nodes: set[str] = set()
        while q:
            side, node = q.popleft()
            if side == "L":
                if node in visited_left:
                    continue
                visited_left.add(node)
                left_nodes.add(node)
                for right_id in sorted(left_to_right[node]):
                    q.append(("R", right_id))
            else:
                if node in visited_right:
                    continue
                visited_right.add(node)
                right_nodes.add(node)
                for left_id in sorted(right_to_left[node]):
                    q.append(("L", left_id))
        component_edges = [
            edge_by_pair[(l_id, r_id)]
            for l_id in sorted(left_nodes)
            for r_id in sorted(right_nodes)
            if (l_id, r_id) in edge_by_pair
        ]
        components.append((sorted(left_nodes), sorted(right_nodes), component_edges))
    return components


def select_one_to_one(edges: list[dict[str, object]]) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    for left_nodes, right_nodes, component_edges in connected_components(edges):
        edge_by_pair = {(e["left_block_id"], e["right_block_id"]): e for e in component_edges}
        if len(left_nodes) * len(right_nodes) <= 400:
            selected_pairs = _exact_component_assignment(left_nodes, right_nodes, edge_by_pair)
        else:
            selected_pairs = _greedy_component_assignment(component_edges)
        for pair in selected_pairs:
            selected.append(edge_by_pair[pair])
    return sorted(selected, key=lambda e: (str(e["left_block_id"]), str(e["right_block_id"])))


def _exact_component_assignment(left_nodes: list[str], right_nodes: list[str],
                                edge_by_pair: dict[tuple[str, str], dict[str, object]]) -> list[tuple[str, str]]:
    if len(left_nodes) <= len(right_nodes):
        rows, cols, transposed = left_nodes, right_nodes, False
    else:
        rows, cols, transposed = right_nodes, left_nodes, True

    @lru_cache(maxsize=None)
    def best(i: int, used_mask: int) -> tuple[tuple[int, float, float], tuple[tuple[str, str], ...]]:
        if i == len(rows):
            return (0, 0.0, 0.0), ()
        skip_score, skip_pairs = best(i + 1, used_mask)
        best_score, best_pairs = skip_score, skip_pairs
        for j, col in enumerate(cols):
            if used_mask & (1 << j):
                continue
            pair = (col, rows[i]) if transposed else (rows[i], col)
            edge = edge_by_pair.get(pair)
            if not edge:
                continue
            tail_score, tail_pairs = best(i + 1, used_mask | (1 << j))
            score = (
                tail_score[0] + int(edge["intersection_bp"]),
                tail_score[1] + float(edge["jaccard"]),
                tail_score[2] + float(edge["left_coverage"]) + float(edge["right_coverage"]),
            )
            pairs = tuple(sorted((pair,) + tail_pairs))
            if score > best_score or (score == best_score and pairs < best_pairs):
                best_score, best_pairs = score, pairs
        return best_score, best_pairs

    return list(best(0, 0)[1])


def _greedy_component_assignment(edges: list[dict[str, object]]) -> list[tuple[str, str]]:
    used_left: set[str] = set()
    used_right: set[str] = set()
    selected = []
    for edge in sorted(edges, key=lambda e: (-int(e["intersection_bp"]), -float(e["jaccard"]),
                                            str(e["left_block_id"]), str(e["right_block_id"]))):
        left_id, right_id = str(edge["left_block_id"]), str(edge["right_block_id"])
        if left_id in used_left or right_id in used_right:
            continue
        used_left.add(left_id)
        used_right.add(right_id)
        selected.append((left_id, right_id))
    return selected


def classify_relationships(left: dict[str, Family], right: dict[str, Family],
                           edges: list[dict[str, object]]) -> list[dict[str, object]]:
    left_adj: dict[str, set[str]] = defaultdict(set)
    right_adj: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        left_adj[str(edge["left_block_id"])].add(str(edge["right_block_id"]))
        right_adj[str(edge["right_block_id"])].add(str(edge["left_block_id"]))
    rows = []
    seen_left: set[str] = set()
    seen_right: set[str] = set()
    for left_nodes, right_nodes, component_edges in connected_components(edges):
        seen_left.update(left_nodes)
        seen_right.update(right_nodes)
        if len(left_nodes) == 1 and len(right_nodes) == 1:
            rel = "one_to_one"
        elif len(left_nodes) == 1:
            rel = "split"
        elif len(right_nodes) == 1:
            rel = "merge"
        else:
            rel = "complex"
        rows.append({
            "relationship": rel,
            "left_block_ids": ",".join(left_nodes),
            "right_block_ids": ",".join(right_nodes),
            "n_left": len(left_nodes),
            "n_right": len(right_nodes),
            "n_edges": len(component_edges),
            "intersection_bp": sum(int(edge["intersection_bp"]) for edge in component_edges),
        })
    for left_id in sorted(set(left) - seen_left):
        rows.append({"relationship": "left_only", "left_block_ids": left_id, "right_block_ids": "",
                     "n_left": 1, "n_right": 0, "n_edges": 0, "intersection_bp": 0})
    for right_id in sorted(set(right) - seen_right):
        rows.append({"relationship": "right_only", "left_block_ids": "", "right_block_ids": right_id,
                     "n_left": 0, "n_right": 1, "n_edges": 0, "intersection_bp": 0})
    order = {"one_to_one": 0, "split": 1, "merge": 2, "complex": 3, "left_only": 4, "right_only": 5}
    return sorted(rows, key=lambda r: (order[str(r["relationship"])], str(r["left_block_ids"]), str(r["right_block_ids"])))


def consensus_intervals(family: Family) -> dict[str, list[tuple[int, int]]]:
    return {chrom: merge_intervals(runs) for chrom, runs in sorted(intervals_by_chrom(family.intervals).items())}


def family_metrics(families: dict[str, Family], consistency_threshold: float) -> list[dict[str, object]]:
    rows = []
    for block_id in sorted(families):
        family = families[block_id]
        consensus = consensus_intervals(family)
        consensus_bp = interval_length(consensus)
        status_counts = Counter(family.occurrence_status.values())
        occurrence_jaccards: list[float] = []
        occurrence_precisions: list[float] = []
        occurrence_recalls: list[float] = []
        orientation_values = [i.strand for i in family.intervals if i.strand in {"+", "-"}]
        dominant = dominant_strand(orientation_values)
        orientation_consistent = (
            sum(1 for s in orientation_values if s == dominant) / len(orientation_values)
            if orientation_values and dominant != "." else None
        )
        by_occ: dict[str, list[Interval]] = defaultdict(list)
        for interval in family.intervals:
            by_occ[interval.occurrence_id].append(interval)
        passing = 0
        for occurrence_id in sorted(family.occurrence_ids):
            occ_intervals = intervals_by_chrom(by_occ.get(occurrence_id, []))
            occ_bp = interval_length(occ_intervals)
            inter = intersection_length(occ_intervals, consensus)
            union = occ_bp + consensus_bp - inter
            precision = inter / occ_bp if occ_bp else 0.0
            recall = inter / consensus_bp if consensus_bp else 0.0
            jaccard = inter / union if union else 1.0
            occurrence_precisions.append(precision)
            occurrence_recalls.append(recall)
            occurrence_jaccards.append(jaccard)
            passing += int(jaccard >= consistency_threshold)
        mapped_occurrences = len([oid for oid in family.occurrence_ids if family.occurrence_status.get(oid) != "unmapped"])
        rows.append({
            "block_id": block_id,
            "n_occurrences": len(family.occurrence_ids),
            "n_species": len(family.species),
            "unique_occurrences": status_counts["unique"],
            "duplicated_occurrences": status_counts["duplicated"],
            "ambiguous_occurrences": status_counts["ambiguous"],
            "unmapped_occurrences": status_counts["unmapped"],
            "fraction_unique": _ratio(status_counts["unique"], len(family.occurrence_ids)),
            "fraction_duplicated": _ratio(status_counts["duplicated"], len(family.occurrence_ids)),
            "fraction_ambiguous": _ratio(status_counts["ambiguous"], len(family.occurrence_ids)),
            "fraction_unmapped": _ratio(status_counts["unmapped"], len(family.occurrence_ids)),
            "ancestral_chromosomes": ",".join(sorted(consensus)),
            "consensus_intervals": _encode_consensus(consensus),
            "consensus_bp": consensus_bp,
            "mean_occurrence_precision": f"{_mean(occurrence_precisions):.6f}",
            "mean_occurrence_recall": f"{_mean(occurrence_recalls):.6f}",
            "mean_occurrence_jaccard": f"{_mean(occurrence_jaccards):.6f}",
            "family_homology_purity": f"{mapped_occurrences / len(family.occurrence_ids):.6f}" if family.occurrence_ids else "0.000000",
            "orientation_consistency": _fmt_optional(orientation_consistent),
            "fraction_passing_consistency": f"{_ratio(passing, len(family.occurrence_ids)):.6f}",
        })
    return rows


def _ratio(num: int, den: int) -> float:
    return num / den if den else 0.0


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _encode_consensus(consensus: dict[str, list[tuple[int, int]]]) -> str:
    parts = []
    for chrom in sorted(consensus):
        parts.extend(f"{chrom}:{start}-{end}" for start, end in consensus[chrom])
    return ",".join(parts)


EDGE_FIELDS = (
    "left_block_id", "right_block_id", "left_bp", "right_bp", "intersection_bp",
    "union_bp", "jaccard", "left_coverage", "right_coverage",
    "orientation_agreement", "supporting_species", "fraction_unique_support",
    "fraction_duplicated_or_ambiguous",
)
MATCH_FIELDS = EDGE_FIELDS
REL_FIELDS = (
    "relationship", "left_block_ids", "right_block_ids", "n_left", "n_right",
    "n_edges", "intersection_bp",
)
FAMILY_FIELDS = (
    "block_id", "n_occurrences", "n_species", "unique_occurrences",
    "duplicated_occurrences", "ambiguous_occurrences", "unmapped_occurrences",
    "fraction_unique", "fraction_duplicated", "fraction_ambiguous",
    "fraction_unmapped", "ancestral_chromosomes", "consensus_intervals",
    "consensus_bp", "mean_occurrence_precision", "mean_occurrence_recall",
    "mean_occurrence_jaccard", "family_homology_purity",
    "orientation_consistency", "fraction_passing_consistency",
)


def match_mapped(left_mapped: Path, right_mapped: Path, output_prefix: Path,
                 thresholds: dict[str, float], consistency_threshold: float = 0.5) -> dict[str, object]:
    left, left_ancestor = load_families(left_mapped)
    right, right_ancestor = load_families(right_mapped)
    if left_ancestor and right_ancestor and left_ancestor != right_ancestor:
        raise ValueError(f"inputs use different ancestors: {left_ancestor!r} vs {right_ancestor!r}")
    ancestor = left_ancestor or right_ancestor
    edges = candidate_edges(left, right, thresholds)
    matches = select_one_to_one(edges)
    relationships = classify_relationships(left, right, edges)
    left_metrics = family_metrics(left, consistency_threshold)
    right_metrics = family_metrics(right, consistency_threshold)
    prefix = output_prefix
    write_tsv(prefix.with_suffix(".overlap_edges.tsv"), EDGE_FIELDS, edges)
    write_tsv(prefix.with_suffix(".block_matches.tsv"), MATCH_FIELDS, matches)
    write_tsv(prefix.with_suffix(".relationships.tsv"), REL_FIELDS, relationships)
    write_tsv(prefix.with_suffix(".left_family_metrics.tsv"), FAMILY_FIELDS, left_metrics)
    write_tsv(prefix.with_suffix(".right_family_metrics.tsv"), FAMILY_FIELDS, right_metrics)
    summary = {
        "package_version": __version__,
        "ancestor": ancestor,
        "left_mapped": str(left_mapped),
        "right_mapped": str(right_mapped),
        "thresholds": thresholds,
        "consistency_threshold": consistency_threshold,
        "counts": {
            "left_families": len(left),
            "right_families": len(right),
            "overlap_edges": len(edges),
            "one_to_one_matches": len(matches),
            "relationships": Counter(row["relationship"] for row in relationships),
        },
    }
    with prefix.with_suffix(".summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="hal-site-test-") as tmp:
        root = Path(tmp)
        pred = root / "pred.tsv"
        truth = root / "truth.tsv"
        out = root / "scores.tsv"
        common = {"block_id": "B1", "species": "eagle", "ancestor": "avianAnc",
                  "anc_chrom": "Anc4", "strand": "-"}
        write_tsv(pred, RUN_FIELDS, [
            {**common, "anc_start": "120", "anc_end": "200"},
            {**common, "anc_start": "250", "anc_end": "280"},
        ])
        write_tsv(truth, RUN_FIELDS, [
            {**common, "anc_start": "100", "anc_end": "200"},
            {**common, "anc_start": "250", "anc_end": "300"},
        ])
        summary = score_runs(pred, truth, out)
        assert abs(float(summary["micro_precision"]) - 1.0) < 1e-12
        assert abs(float(summary["micro_recall"]) - 110 / 150) < 1e-12
        assert abs(float(summary["micro_jaccard"]) - 110 / 150) < 1e-12
        assert summary["orientation_accuracy"] == 1.0

        left = root / "left.tsv"
        right = root / "right.tsv"
        write_tsv(left, MAPPED_FIELDS, [
            {"block_id": "L1", "occurrence_id": "lo1", "species": "eagle", "ancestor": "avianAnc",
             "anc_chrom": "Anc4", "anc_start": 100, "anc_end": 200, "strand": "+",
             "copy_id": "1", "mapping_status": "unique"},
        ])
        write_tsv(right, MAPPED_FIELDS, [
            {"block_id": "R9", "occurrence_id": "ro1", "species": "eagle", "ancestor": "avianAnc",
             "anc_chrom": "Anc4", "anc_start": 100, "anc_end": 200, "strand": "+",
             "copy_id": "1", "mapping_status": "unique"},
        ])
        summary = match_mapped(left, right, root / "match", default_thresholds())
        assert summary["counts"]["overlap_edges"] == 1
    print("SELF-TEST PASSED")


def default_thresholds() -> dict[str, float]:
    return {
        "min_intersection_bp": 1,
        "min_jaccard": 0.0,
        "min_left_coverage": 0.0,
        "min_right_coverage": 0.0,
        "min_supporting_species": 0,
    }


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    m = sub.add_parser("map", help="map block intervals to one HAL ancestral genome")
    m.add_argument("--input", type=Path, required=True)
    m.add_argument("--hal", type=Path, required=True)
    m.add_argument("--ancestor", required=True)
    m.add_argument("--output", type=Path, required=True)
    m.add_argument("--hal-liftover", default="halLiftover")
    s = sub.add_parser("score", help="legacy same-label comparison of two mapped run tables")
    s.add_argument("--pred", type=Path, required=True)
    s.add_argument("--truth", type=Path, required=True)
    s.add_argument("--output", type=Path, required=True)
    mt = sub.add_parser("match", help="label-independent comparison of two mapped block systems")
    mt.add_argument("--left-mapped", type=Path, required=True)
    mt.add_argument("--right-mapped", type=Path, required=True)
    mt.add_argument("--output-prefix", type=Path, required=True)
    mt.add_argument("--min-intersection-bp", type=int, default=1)
    mt.add_argument("--min-jaccard", type=float, default=0.0)
    mt.add_argument("--min-left-coverage", type=float, default=0.0)
    mt.add_argument("--min-right-coverage", type=float, default=0.0)
    mt.add_argument("--min-supporting-species", type=int, default=0)
    mt.add_argument("--consistency-threshold", type=float, default=0.5)
    sub.add_parser("self-test", help="test interval scoring and matching without HAL")
    return p


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "map":
            map_to_ancestor(args.input, args.hal, args.ancestor, args.output, args.hal_liftover)
        elif args.command == "score":
            print(json.dumps(score_runs(args.pred, args.truth, args.output), indent=2, sort_keys=True))
        elif args.command == "match":
            thresholds = {
                "min_intersection_bp": args.min_intersection_bp,
                "min_jaccard": args.min_jaccard,
                "min_left_coverage": args.min_left_coverage,
                "min_right_coverage": args.min_right_coverage,
                "min_supporting_species": args.min_supporting_species,
            }
            print(json.dumps(match_mapped(
                args.left_mapped, args.right_mapped, args.output_prefix,
                thresholds, args.consistency_threshold,
            ), indent=2, sort_keys=True))
        elif args.command == "self-test":
            self_test()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
