# HAL Site Evaluator

`hal-site-evaluator` evaluates predicted cross-species synteny blocks against a
Cactus/HAL alignment without assigning and storing a separate integer ID at
every nucleotide.

It maps predicted and reference intervals to a common ancestral genome, keeps
the mappings as run-length intervals, and computes overlap using interval
arithmetic. This makes the evaluation scale with the number of mapped alignment
runs rather than the total number of bases.

## What it measures

For every `(species, block_id)` occurrence, the evaluator reports:

- homology precision;
- homology recall;
- homology Jaccard similarity;
- orientation agreement;
- predicted, reference, intersection, and union lengths in ancestral bases.

The current release assumes that predictions and reference data use the same
block-ID vocabulary. It evaluates block recovery; it does not yet match
arbitrarily permuted predicted labels to reference families.

## Requirements

- Python 3.9 or newer
- A HAL alignment
- `halLiftover` on `PATH`

`halLiftover` and `halStats` are command-line programs distributed with the
[Comparative Genomics Toolkit HAL package](https://github.com/ComparativeGenomicsToolkit/hal).
They are external executables and are not installed by this Python package.

## Installation

Clone the repository and install it into a virtual environment:

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

## Input format

Prediction and truth files are tab-separated with a header:

```text
block_id\tspecies\tchrom\tstart\tend\tstrand
B17\tchicken\tchr4\t1000\t2500\t+
B17\teagle\tchr4\t4300\t5790\t-
```

Coordinates must be zero-based and half-open, as in BED. Within each input
file, `(species, block_id)` must be unique. Use a unique block ID for each
separate occurrence if a genome contains multiple copies.

## Choose the ancestral coordinate system

Inspect a HAL file with:

```bash
halStats alignment.hal
halStats alignment.hal --genomes
halStats alignment.hal --tree
```

For a bird-only experiment, choose the avian ancestral genome represented in
the HAL tree. Use its exact HAL genome name as `--ancestor`.

## Usage

Map predictions to the chosen ancestor:

```bash
hal-site-eval map \
  --input predictions.tsv \
  --hal vgp.hal \
  --ancestor avianAncestor \
  --output predictions.ancestral.tsv
```

Map reference blocks in the same way:

```bash
hal-site-eval map \
  --input truth.tsv \
  --hal vgp.hal \
  --ancestor avianAncestor \
  --output truth.ancestral.tsv
```

Score the mapped runs:

```bash
hal-site-eval score \
  --pred predictions.ancestral.tsv \
  --truth truth.ancestral.tsv \
  --output block_scores.tsv
```

The aggregate summary is printed as JSON. Per-occurrence results are written
to `block_scores.tsv`.

If `halLiftover` is installed under another name or path, specify it explicitly:

```bash
hal-site-eval map ... --hal-liftover /path/to/halLiftover
```

## Testing

### Test interval scoring without HAL

```bash
hal-site-eval self-test
pytest -q
```

### Smoke-test a real HAL file

1. Select one genome and 10--20 short, uniquely aligned intervals.
2. Use the same TSV as both prediction and truth.
3. Map both files to the same ancestor and score them.
4. Precision, recall, Jaccard, and orientation accuracy should all be `1.0`.
5. Shift a boundary, remove an interval, add an unrelated interval, and reverse
   a strand. Confirm that recall, precision, and orientation respond as expected.

### Measure scaling

```bash
/usr/bin/time -v hal-site-eval map \
  --input predictions.tsv \
  --hal vgp.hal \
  --ancestor avianAncestor \
  --output predictions.ancestral.tsv
```

Repeat with 1, 5, 10, and 20 species and record elapsed time and peak resident
memory. The mapper invokes `halLiftover` once per species, not once per block.

## Interpretation and current limitations

- HAL is treated as evaluation truth; the package does not assess whether the
  underlying Cactus alignment is biologically correct.
- HAL provides nucleotide homology but does not define a unique biological
  synteny-block segmentation.
- Unmapped sequence contributes zero overlap.
- Multiple lifted copies are currently unioned. A future version should report
  unique, duplicated, ambiguous, and unmapped mappings separately.
- Exact block-ID matching is required in version 0.1.0.

## Build a distribution

```bash
python -m build
python -m pip install dist/hal_site_evaluator-0.1.0-py3-none-any.whl
hal-site-eval self-test
```

## License

MIT
