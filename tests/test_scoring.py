from pathlib import Path

import pytest

from hal_site_evaluator import (
    MAPPED_FIELDS,
    RUN_FIELDS,
    SOURCE_FIELDS,
    default_thresholds,
    intersection_length,
    load_families,
    map_to_ancestor,
    match_mapped,
    merge_intervals,
    score_runs,
    validate_source,
    write_tsv,
)


def test_merge_intervals_merges_overlaps_and_touching_runs():
    assert merge_intervals([(10, 20), (15, 25), (25, 30), (40, 50)]) == [
        (10, 30),
        (40, 50),
    ]


def test_intersection_across_chromosomes():
    left = {"Anc1": [(0, 10)], "Anc2": [(0, 20)]}
    right = {"Anc1": [(5, 15)], "Anc2": [(10, 30)]}
    assert intersection_length(left, right) == 15


def test_score_runs(tmp_path: Path):
    pred = tmp_path / "pred.tsv"
    truth = tmp_path / "truth.tsv"
    scores = tmp_path / "scores.tsv"

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

    result = score_runs(pred, truth, scores)
    assert result["micro_precision"] == 1.0
    assert abs(result["micro_recall"] - 110 / 150) < 1e-12
    assert abs(result["micro_jaccard"] - 110 / 150) < 1e-12
    assert result["orientation_accuracy"] == 1.0


def mapped_row(block, occ, species, chrom, start, end, strand="+", status="unique", ancestor="anc"):
    return {
        "block_id": block,
        "occurrence_id": occ,
        "species": species,
        "ancestor": ancestor,
        "anc_chrom": chrom,
        "anc_start": start,
        "anc_end": end,
        "strand": strand,
        "copy_id": "1",
        "mapping_status": status,
        "source_status": "unspecified",
    }


def read_dicts(path: Path):
    import csv

    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_match(tmp_path: Path, left_rows, right_rows, **thresholds):
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    write_tsv(left, MAPPED_FIELDS, left_rows)
    write_tsv(right, MAPPED_FIELDS, right_rows)
    params = default_thresholds()
    params.update(thresholds)
    summary = match_mapped(left, right, tmp_path / "comparison", params)
    return summary, {
        "edges": read_dicts(tmp_path / "comparison.overlap_edges.tsv"),
        "matches": read_dicts(tmp_path / "comparison.block_matches.tsv"),
        "relationships": read_dicts(tmp_path / "comparison.relationships.tsv"),
        "left_metrics": read_dicts(tmp_path / "comparison.left_family_metrics.tsv"),
        "right_metrics": read_dicts(tmp_path / "comparison.right_family_metrics.tsv"),
    }


def test_extended_source_defaults_allow_duplicate_species_block():
    rows = [
        {"block_id": "B", "species": "sp1", "chrom": "chr1", "start": "0", "end": "10", "strand": "+"},
        {"block_id": "B", "species": "sp1", "chrom": "chr1", "start": "20", "end": "30", "strand": "+"},
    ]
    occurrences = validate_source(rows)
    assert len({o.occurrence_id for o in occurrences}) == 2
    assert [o.copy_id for o in occurrences] == ["1", "1"]
    assert [o.status for o in occurrences] == ["unspecified", "unspecified"]


def test_identical_regions_with_different_labels_match(tmp_path: Path):
    summary, outputs = run_match(
        tmp_path,
        [mapped_row("maf_A", "l1", "sp1", "anc1", 0, 100)],
        [mapped_row("hal_99", "r1", "sp1", "anc1", 0, 100)],
    )
    assert summary["counts"]["overlap_edges"] == 1
    assert outputs["edges"][0]["left_block_id"] == "maf_A"
    assert outputs["edges"][0]["right_block_id"] == "hal_99"
    assert outputs["edges"][0]["jaccard"] == "1.000000"
    assert outputs["relationships"][0]["relationship"] == "one_to_one"


def test_partial_interval_overlap_metrics(tmp_path: Path):
    _, outputs = run_match(
        tmp_path,
        [mapped_row("L", "l1", "sp1", "anc1", 0, 100)],
        [mapped_row("R", "r1", "sp1", "anc1", 50, 150)],
    )
    edge = outputs["edges"][0]
    assert edge["intersection_bp"] == "50"
    assert edge["union_bp"] == "150"
    assert edge["left_coverage"] == "0.500000"
    assert edge["right_coverage"] == "0.500000"


def test_split_merge_complex_and_only_classifications(tmp_path: Path):
    left_rows = [
        mapped_row("Lsplit", "ls", "sp1", "anc1", 0, 100),
        mapped_row("Lmerge1", "lm1", "sp1", "anc2", 0, 50),
        mapped_row("Lmerge2", "lm2", "sp1", "anc2", 50, 100),
        mapped_row("Lc1", "lc1", "sp1", "anc3", 0, 100),
        mapped_row("Lc2", "lc2", "sp1", "anc3", 50, 150),
        mapped_row("Lonly", "lo", "sp1", "anc9", 0, 10),
    ]
    right_rows = [
        mapped_row("Rsplit1", "rs1", "sp1", "anc1", 0, 50),
        mapped_row("Rsplit2", "rs2", "sp1", "anc1", 50, 100),
        mapped_row("Rmerge", "rm", "sp1", "anc2", 0, 100),
        mapped_row("Rc1", "rc1", "sp1", "anc3", 25, 75),
        mapped_row("Rc2", "rc2", "sp1", "anc3", 75, 125),
        mapped_row("Ronly", "ro", "sp1", "anc8", 0, 10),
    ]
    _, outputs = run_match(tmp_path, left_rows, right_rows)
    rels = {row["relationship"] for row in outputs["relationships"]}
    assert {"split", "merge", "complex", "left_only", "right_only"} <= rels


def test_duplicate_ambiguous_unmapped_family_metrics(tmp_path: Path):
    _, outputs = run_match(
        tmp_path,
        [
            mapped_row("L", "l_unique", "sp1", "anc1", 0, 100, status="unique"),
            mapped_row("L", "l_dup", "sp2", "anc1", 0, 100, status="duplicated"),
            mapped_row("L", "l_amb", "sp3", "anc2", 0, 50, status="ambiguous"),
            {**mapped_row("L", "l_unmap", "sp4", "", "", "", status="unmapped"), "anc_chrom": ""},
        ],
        [mapped_row("R", "r1", "sp1", "anc1", 0, 100)],
    )
    metrics = outputs["left_metrics"][0]
    assert metrics["n_occurrences"] == "4"
    assert metrics["unique_occurrences"] == "1"
    assert metrics["duplicated_occurrences"] == "1"
    assert metrics["ambiguous_occurrences"] == "1"
    assert metrics["unmapped_occurrences"] == "1"


def test_orientation_agreement_and_mismatch(tmp_path: Path):
    _, outputs = run_match(
        tmp_path,
        [
            mapped_row("Lsame", "ls", "sp1", "anc1", 0, 50, "+"),
            mapped_row("Ldiff", "ld", "sp1", "anc1", 100, 150, "+"),
        ],
        [
            mapped_row("Rsame", "rs", "sp1", "anc1", 0, 50, "+"),
            mapped_row("Rdiff", "rd", "sp1", "anc1", 100, 150, "-"),
        ],
    )
    by_pair = {(e["left_block_id"], e["right_block_id"]): e for e in outputs["edges"]}
    assert by_pair[("Lsame", "Rsame")]["orientation_agreement"] == "1.000000"
    assert by_pair[("Ldiff", "Rdiff")]["orientation_agreement"] == "0.000000"


def test_different_ancestral_chromosomes_do_not_match(tmp_path: Path):
    summary, outputs = run_match(
        tmp_path,
        [mapped_row("L", "l1", "sp1", "anc1", 0, 100)],
        [mapped_row("R", "r1", "sp1", "anc2", 0, 100)],
    )
    assert summary["counts"]["overlap_edges"] == 0
    assert outputs["relationships"][0]["relationship"] == "left_only"
    assert outputs["relationships"][1]["relationship"] == "right_only"


def test_multi_run_interval_family_consensus(tmp_path: Path):
    _, outputs = run_match(
        tmp_path,
        [
            mapped_row("L", "l1", "sp1", "anc1", 0, 10),
            mapped_row("L", "l1", "sp1", "anc1", 20, 30),
        ],
        [mapped_row("R", "r1", "sp1", "anc1", 0, 30)],
    )
    assert outputs["left_metrics"][0]["consensus_intervals"] == "anc1:0-10,anc1:20-30"


def test_deterministic_tied_match_resolution(tmp_path: Path):
    _, outputs = run_match(
        tmp_path,
        [mapped_row("L1", "l1", "sp1", "anc1", 0, 100), mapped_row("L2", "l2", "sp1", "anc1", 0, 100)],
        [mapped_row("R1", "r1", "sp1", "anc1", 0, 100), mapped_row("R2", "r2", "sp1", "anc1", 0, 100)],
    )
    assert [(m["left_block_id"], m["right_block_id"]) for m in outputs["matches"]] == [
        ("L1", "R1"),
        ("L2", "R2"),
    ]


def test_threshold_filtering(tmp_path: Path):
    summary, outputs = run_match(
        tmp_path,
        [mapped_row("L", "l1", "sp1", "anc1", 0, 100)],
        [mapped_row("R", "r1", "sp1", "anc1", 50, 150)],
        min_jaccard=0.5,
    )
    assert summary["counts"]["overlap_edges"] == 0
    assert outputs["matches"] == []


def test_empty_inputs(tmp_path: Path):
    summary, outputs = run_match(tmp_path, [], [])
    assert summary["counts"]["left_families"] == 0
    assert outputs["edges"] == []


def test_malformed_coordinates_rejected(tmp_path: Path):
    left = tmp_path / "bad.tsv"
    right = tmp_path / "right.tsv"
    write_tsv(left, MAPPED_FIELDS, [mapped_row("L", "l1", "sp1", "anc1", 10, 0)])
    write_tsv(right, MAPPED_FIELDS, [])
    with pytest.raises(ValueError):
        match_mapped(left, right, tmp_path / "out", default_thresholds())


def test_legacy_mapped_format_still_loads(tmp_path: Path):
    legacy = tmp_path / "legacy.tsv"
    write_tsv(legacy, RUN_FIELDS, [
        {"block_id": "B", "species": "sp1", "ancestor": "anc", "anc_chrom": "anc1",
         "anc_start": "0", "anc_end": "10", "strand": "+"}
    ])
    families, ancestor = load_families(legacy)
    assert ancestor == "anc"
    assert families["B"].occurrence_status


def test_mocked_liftover_preserves_occurrence_metadata(tmp_path: Path, monkeypatch):
    source = tmp_path / "source.tsv"
    output = tmp_path / "mapped.tsv"
    write_tsv(source, SOURCE_FIELDS, [
        {"block_id": "B", "occurrence_id": "occ1", "species": "sp1", "chrom": "chr1",
         "start": "0", "end": "100", "strand": "+", "copy_id": "2", "status": "candidate"}
    ])
    monkeypatch.setattr("hal_site_evaluator.shutil.which", lambda exe: exe)

    def fake_run(cmd, text, capture_output):
        target = Path(cmd[-1])
        with target.open("w") as handle:
            handle.write("anc1\t10\t110\tHSE000000000001\t0\t+\n")

        class Completed:
            returncode = 0
            stderr = ""

        return Completed()

    monkeypatch.setattr("hal_site_evaluator.subprocess.run", fake_run)
    map_to_ancestor(source, tmp_path / "fake.hal", "anc", output, "halLiftover")
    rows = read_dicts(output)
    assert rows[0]["block_id"] == "B"
    assert rows[0]["occurrence_id"] == "occ1"
    assert rows[0]["copy_id"] == "2"
    assert rows[0]["source_status"] == "candidate"
    assert rows[0]["mapping_status"] == "unique"
