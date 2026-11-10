# Example

This directory contains a small worked example for `hal-site-evaluator` using two
synthetic block sets (`left_blocks` and `right_blocks`) with two species (human, mouse)
and a toy ancestral genome called `Ancestor`.

## Files

| File | Description |
|---|---|
| `left_blocks.tsv` | Input block occurrence table — 3 blocks (L1, L2, L3) |
| `right_blocks.tsv` | Input block occurrence table — 4 blocks (R1, R2, R3, R4) |
| `left_blocks.mapped.tsv` | Output of `map` on `left_blocks.tsv` |
| `right_blocks.mapped.tsv` | Output of `map` on `right_blocks.tsv` |
| `results.*` | Output of `match` on the two mapped files |

The designed relationships between the two block sets are:

| Relationship | Left | Right |
|---|---|---|
| one_to_one | L1 | R1 |
| split | L2 | R2, R3 |
| left_only | L3 | — |
| right_only | — | R4 |

## Commands

The mapped TSVs and results were generated with:

```bash
# Step 1 — map each block set to the ancestral genome
hal-site-eval map \
  --input left_blocks.tsv \
  --hal alignment.hal \
  --ancestor Ancestor \
  --output left_blocks.mapped.tsv

hal-site-eval map \
  --input right_blocks.tsv \
  --hal alignment.hal \
  --ancestor Ancestor \
  --output right_blocks.mapped.tsv

# Step 2 — compare the two mapped sets
hal-site-eval match \
  --left-mapped left_blocks.mapped.tsv \
  --right-mapped right_blocks.mapped.tsv \
  --output-prefix results
```

Or equivalently, using the `eval` command which does both steps in one go:

```bash
hal-site-eval eval \
  --left left_blocks.tsv \
  --right right_blocks.tsv \
  --hal alignment.hal \
  --ancestor Ancestor \
  --output-prefix results \
  --keep-mapped
```

> **Note:** The mapped TSVs here were pre-computed with simulated ancestral coordinates
> rather than a real HAL file. To run the commands above yourself you need an actual
> `alignment.hal` file with a genome named `Ancestor` that contains human and mouse.
