from pathlib import Path

from hal_site_evaluator import (
    RUN_FIELDS,
    intersection_length,
    merge_intervals,
    score_runs,
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

