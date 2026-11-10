# HAL Site Evaluator

**HAL Site Evaluator** (`hal-site-evaluator`) is a tool for evaluating and comparing cross-species synteny block systems. It maps block occurrence intervals into a shared ancestral coordinate system defined by a [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) HAL alignment, then compares the resulting intervals to assess consistency and agreement between methods.

## Requirements

- Python 3.9+
- [`halLiftover`](https://github.com/ComparativeGenomicsToolkit/hal) on `PATH`

## Installation

```bash
git clone https://github.com/ytabatabaee/hal-site-evaluator.git
cd hal-site-evaluator
pip install -e .
```

Verify the installation:

```bash
hal-site-eval self-test
```

## Usage

### Quick start

Pass two block sets directly and get results in one command:

```bash
hal-site-eval eval \
  --left blocks_A.tsv \
  --right blocks_B.tsv \
  --hal alignment.hal \
  --ancestor Ancestor \
  --output-prefix results/comparison
```

This maps both inputs to the ancestral genome internally, then runs the comparison. Add `--keep-mapped` to also save the intermediate mapped TSVs (as `results/comparison.left_mapped.tsv` and `results/comparison.right_mapped.tsv`).

Optional filters:

| Option | Description |
|---|---|
| `--min-intersection-bp` | Minimum overlapping base pairs |
| `--min-jaccard` | Minimum Jaccard score |
| `--min-left-coverage` | Minimum coverage on left side |
| `--min-right-coverage` | Minimum coverage on right side |
| `--min-supporting-species` | Minimum number of supporting species |

### Advanced: map and match separately

If you are comparing multiple block systems against each other, map each one once and reuse the mapped file across comparisons:

```bash
hal-site-eval map --input blocks_A.tsv --hal alignment.hal --ancestor Ancestor --output blocks_A.mapped.tsv
hal-site-eval map --input blocks_B.tsv --hal alignment.hal --ancestor Ancestor --output blocks_B.mapped.tsv

hal-site-eval match \
  --left-mapped blocks_A.mapped.tsv \
  --right-mapped blocks_B.mapped.tsv \
  --output-prefix results/comparison
```

`match` accepts the same optional filters as `eval`.

## Input Format

Tab-separated occurrence table, one row per block occurrence:

```
block_id    species    chrom    start    end    strand
```

Coordinates are **zero-based, half-open**. Three optional columns can be appended: `occurrence_id`, `copy_id`, `status`.

Example:

```
block_id	species	chrom	start	end	strand
L1	human	chr1	1000	2000	+
L1	mouse	chr1	1100	2100	+
L2	human	chr3	5000	7000	-
L2	mouse	chr3	4800	6800	-
L3	human	chr5	100	500	+
```

Each row is one occurrence of the block in one species. See [`example/`](example/) for a complete worked example.

## Output Files

`map` produces a mapped TSV with an added `anc_intervals` column and a per-occurrence status (`unique`, `duplicated`, `ambiguous`, or `unmapped`).

`match` produces the following files under the given prefix:

| File | Contents |
|---|---|
| `PREFIX.overlap_edges.tsv` | Pairwise overlaps between left and right blocks |
| `PREFIX.block_matches.tsv` | One-to-one best matches |
| `PREFIX.relationships.tsv` | Relationship class per connected component (`one_to_one`, `split`, `merge`, `complex`, `left_only`, `right_only`) |
| `PREFIX.left_family_metrics.tsv` | Per-family metrics for the left block set |
| `PREFIX.right_family_metrics.tsv` | Per-family metrics for the right block set |
| `PREFIX.summary.json` | Run parameters and aggregate counts |

## License

MIT
