# HAL Site Evaluator

`hal-site-evaluator` evaluates proposed cross-species synteny block systems by
mapping their intervals into a common Cactus/HAL ancestral coordinate system and
comparing the resulting run-length intervals. It does not assign or store a
separate integer ID at every nucleotide.

HAL is used as the nucleotide-homology reference for these evaluations. A block
system such as MAF2Synteny output is treated as one method for segmenting that
homology into block families, not as biological ground truth. The package can
therefore compare MAF2Synteny blocks, HAL-derived candidate blocks, or any
future block extractor that emits the supported occurrence table.

All coordinates are zero-based, half-open intervals.

## What It Measures

The evaluator supports two complementary tasks:

- Intrinsic family consistency: whether occurrences assigned to the same block
  family map to consistent ancestral homologous regions.
- Label-independent method comparison: whether two independently named block
  systems describe the same ancestral homology, including split, merge, and
  complex relationships.

The legacy `score` command is still available for same-label mapped tables. It
reports ancestral precision, recall, Jaccard, orientation agreement, and mapped
interval lengths for matching `(species, block_id)` labels.

## Requirements

- Python 3.9 or newer
- A HAL alignment for `map`
- `halLiftover` on `PATH` for `map`

`halLiftover` and `halStats` are command-line programs distributed with the
[Comparative Genomics Toolkit HAL package](https://github.com/ComparativeGenomicsToolkit/hal).
They are external executables and are not installed by this Python package.

## Installation

```bash
git clone https://github.com/ytabatabaee/hal-site-evaluator.git
cd hal-site-evaluator
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

For development and testing:

```bash
python -m pip install pytest build
```

## Source Occurrence Format

Input to `hal-site-eval map` is tab-separated. The required legacy columns are:

```text
block_id species chrom start end strand
```

Version 0.2.0 also accepts copy-aware occurrence metadata:

```text
block_id occurrence_id species chrom start end strand copy_id status
```

`occurrence_id`, `copy_id`, and `status` are optional for backward
compatibility. When they are missing, the mapper generates a stable deterministic
`occurrence_id`, sets `copy_id` to `1`, and sets `status` to `unspecified`.
More than one occurrence of the same block in the same species is allowed;
occurrence identity, not `(species, block_id)`, is the uniqueness constraint.

Example:

```text
block_id	occurrence_id	species	chrom	start	end	strand	copy_id	status
B17	B17.chicken.1	chicken	chr4	1000	2500	+	1	candidate
B17	B17.chicken.2	chicken	chr9	700	1300	-	2	candidate
B17	B17.eagle.1	eagle	chr4	4300	5790	-	1	candidate
```

## Mapping Semantics

`map` batches records by species and invokes `halLiftover` once per species.
Temporary BED names are internal safe record IDs; block IDs are not used as the
BED record name. The output preserves `block_id`, `occurrence_id`, `copy_id`,
source status, source coordinates, and mapping coverage statistics.

Mapping status is assigned per source occurrence:

- `unique`: all lifted runs are on one ancestral chromosome, have a consistent
  orientation, do not overlap as duplicated lifted coverage, and cover the
  source length.
- `duplicated`: lifted coverage is larger than the source interval or lifted
  runs overlap in ancestral space, indicating multiple candidate copies.
- `ambiguous`: lifted runs span multiple ancestral chromosomes, disagree in
  orientation, or only partially cover the source interval.
- `unmapped`: no lifted ancestral interval was returned.

Split alignment runs that remain on the same ancestral chromosome with
consistent orientation are retained as separate run-length intervals for one
occurrence. HAL BED output alone does not always distinguish paralogous copies
from fragmented alignment confidently; uncertain cases are marked `ambiguous`
instead of being silently collapsed.

Coverage columns report how much source sequence was uniquely mapped,
duplicated, ambiguous, or unmapped. Multiple lifted runs are retained; the
mapper does not select one paralog and does not union biologically distinct
copies into one interval.

## Example Workflow

Both inputs must be mapped to the same ancestral HAL genome:

```bash
hal-site-eval map \
  --input maf2synteny_blocks.tsv \
  --hal vgp.hal \
  --ancestor avianAncestor \
  --output maf2synteny.ancestral.tsv

hal-site-eval map \
  --input hal_extracted_blocks.tsv \
  --hal vgp.hal \
  --ancestor avianAncestor \
  --output hal_blocks.ancestral.tsv

hal-site-eval match \
  --left-mapped maf2synteny.ancestral.tsv \
  --right-mapped hal_blocks.ancestral.tsv \
  --output-prefix maf_vs_hal
```

The terms `left` and `right` are neutral. The command does not require matching
block IDs and does not treat either side as unquestionable truth.

## Match Outputs

`match` writes deterministic, machine-readable files:

- `PREFIX.overlap_edges.tsv`: all thresholded left-right family overlaps.
- `PREFIX.block_matches.tsv`: a deterministic one-to-one summary matching.
- `PREFIX.relationships.tsv`: `one_to_one`, `split`, `merge`, `complex`,
  `left_only`, and `right_only` classifications.
- `PREFIX.left_family_metrics.tsv`: intrinsic metrics for left families.
- `PREFIX.right_family_metrics.tsv`: intrinsic metrics for right families.
- `PREFIX.summary.json`: parameters, thresholds, package version, ancestor, and
  aggregate counts.

Candidate family pairs are generated only from overlapping ancestral intervals
on the same chromosome. Individual bases are never materialized.

Filtering options:

```text
--min-intersection-bp
--min-jaccard
--min-left-coverage
--min-right-coverage
--min-supporting-species
```

## Metric Definitions

Overlap edge metrics:

- `intersection_bp`: ancestral bases shared by a left and right family.
- `union_bp`: left family span plus right family span minus intersection.
- `jaccard`: `intersection_bp / union_bp`.
- `left_coverage`: fraction of left ancestral span overlapped by the right
  family.
- `right_coverage`: fraction of right ancestral span overlapped by the left
  family.
- `orientation_agreement`: overlap-weighted fraction of oriented bases with the
  same strand, or `NA` when orientation is unavailable.
- `supporting_species`: number of species with overlapping left and right runs.
- `fraction_unique_support`: fraction of supporting mapped runs marked
  `unique`.
- `fraction_duplicated_or_ambiguous`: fraction of supporting mapped runs marked
  `duplicated` or `ambiguous`.

Intrinsic family metrics:

- `n_occurrences`: number of occurrence IDs in the block family.
- `n_species`: number of represented species.
- status counts and fractions: unique, duplicated, ambiguous, and unmapped
  occurrences.
- `ancestral_chromosomes`: ancestral chromosomes represented by mapped runs.
- `consensus_intervals`: the family consensus represented as one or more
  run-length intervals, not necessarily one continuous span.
- `consensus_bp`: total consensus interval length.
- `mean_occurrence_precision`: mean fraction of each occurrence covered by the
  family consensus.
- `mean_occurrence_recall`: mean fraction of the consensus covered by each
  occurrence.
- `mean_occurrence_jaccard`: mean occurrence-to-consensus Jaccard.
- `family_homology_purity`: fraction of occurrences with mapped, non-unmapped
  ancestral evidence.
- `orientation_consistency`: fraction of oriented mapped runs matching the
  dominant family orientation, or `NA` when orientation is unavailable.
- `fraction_passing_consistency`: fraction of occurrences whose
  occurrence-to-consensus Jaccard is at least `--consistency-threshold`.

The consensus is the union of mapped family runs in ancestral coordinates. This
is intentionally represented as multiple intervals when the family maps to a
multi-run ancestral region.

## One-To-One Matching and Segmentation

The overlap graph is the authoritative comparison. The one-to-one match table is
only a convenient summary. For each connected component of the overlap graph,
the evaluator chooses a deterministic maximum-weight bipartite assignment. The
primary weight is homologous intersection length; Jaccard and coverage are used
as secondary weights, and label order breaks exact ties deterministically.

Segmentation classifications are component-level:

- `one_to_one`: one left family overlaps one right family.
- `split`: one left family overlaps multiple right families.
- `merge`: multiple left families overlap one right family.
- `complex`: multiple left families overlap multiple right families.
- `left_only`: a left family has no thresholded right overlap.
- `right_only`: a right family has no thresholded left overlap.

## Legacy Same-Label Scoring

Mapped TSV files with the original columns remain valid:

```text
block_id species ancestor anc_chrom anc_start anc_end strand
```

Use:

```bash
hal-site-eval score \
  --pred predictions.ancestral.tsv \
  --truth truth.ancestral.tsv \
  --output block_scores.tsv
```

The aggregate summary is printed as JSON. Per-occurrence results are written to
the requested TSV.

## Scalability

The mapper invokes `halLiftover` once per species, not once per block. Matching
uses sorted interval sweeps to generate candidate family pairs only when
ancestral chromosome intervals overlap, then processes connected overlap
components independently.

Approximate complexity:

- mapping: dominated by HAL liftover per species plus linear parsing of lifted
  BED runs;
- candidate generation: `O(R log R + K)` per ancestral chromosome, where `R`
  is the number of merged family runs and `K` is the number of overlapping
  candidate pairs;
- edge scoring: proportional to interval runs in candidate pairs;
- one-to-one assignment: exact dynamic programming for moderate connected
  components, with deterministic greedy fallback for very large components.

The implementation avoids all-base site-name materialization and avoids
comparing every left family with every right family.

## Limitations

- HAL is treated as the homology reference; this package does not assess
  whether the underlying Cactus alignment is biologically correct.
- HAL homology does not define a unique synteny-block segmentation.
- Mapping status classification is conservative because generic HAL liftover
  BED output may not fully distinguish true paralogous copies from fragmented or
  ambiguous alignment output.
- This repository evaluates proposed block systems. It does not implement a
  multi-genome HAL block-construction algorithm, species-tree inference, or
  MLGO execution. Those belong in a separate extractor such as the proposed
  `hal-multisynteny` package.

## Testing

```bash
hal-site-eval self-test
python -m pytest -q
```

Smoke-test a real HAL file by mapping a small interval table twice to the same
ancestor and confirming that exact self-comparison gives perfect overlap.

## Build a Distribution

```bash
python -m build
python -m pip install dist/hal_site_evaluator-0.2.0-py3-none-any.whl
hal-site-eval self-test
```

## Changelog

### 0.2.0

- Added copy-aware occurrence metadata and mapping status columns.
- Added label-independent `match` command with overlap graph, one-to-one
  summary matches, relationship classifications, family metrics, and JSON
  summary output.
- Preserved legacy `score`, mapped TSV loading, and `self-test` behavior.

### 0.1.0

- Initial run-length HAL ancestral-coordinate mapping and same-label scoring.

## License

MIT
