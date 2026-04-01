# AGENTS.md — Data Architecture & Column Usage Rules

This file describes the parquet dataset design for this project and the strict rules for which columns to use in each analysis context. **Do not deviate from these rules.**

---

## Dataset Overview

Four parquet files are loaded, identified by their `source` column:

| source             | Description                                      | AF origin      |
|--------------------|--------------------------------------------------|----------------|
| `background`       | Background genes, gnomAD allele frequencies      | Real (gnomAD)  |
| `clingen`          | ClinGen haploinsufficient genes, gnomAD AFs      | Real (gnomAD)  |
| `background_null`  | Same background genes, **permuted** AFs          | Permuted/synthetic |
| `clingen_null`     | Same ClinGen genes, **permuted** AFs             | Permuted/synthetic |

---

## Column Rules — CRITICAL

### Real datasets (`background`, `clingen`)

- **Use `vg_predicted`** — total Vg computed from real gnomAD allele frequencies.
- **Use `vg_common`, `vg_rare`** — regional AF-bin breakdowns using real AFs.
- **Use `vg_{region}`** (e.g. `vg_promoter_core`, `vg_distal_upstream`) — per-region Vg using real AFs.
- **NEVER use `_perm` columns from these files for biological comparisons.** The `vg_predicted_perm`, `vg_{region}_perm` columns exist in these files as an artifact of the pipeline but they should NOT be used as the null/permuted expectation in any analysis. They are not the intended null model.

### Synthetic/null datasets (`background_null`, `clingen_null`)

- **Use `vg_predicted_perm`** — total Vg computed from permuted AFs. This is the null expectation for `vg_predicted`.
- **Use `vg_{region}_perm`** (e.g. `vg_promoter_core_perm`, `vg_distal_upstream_perm`) — per-region Vg under permuted AFs. These are the null expectations for per-region observed Vg.
- The `vg_{region}` columns (without `_perm`) in null files are **zeroed out** by design — do not use them.

---

## Correct Patterns for Ratio / Depletion Analyses

### Total Vg depletion (obs / perm)

```
obs  = all_genes where source in ['background', 'clingen']     → column: vg_predicted
perm = all_genes where source in ['background_null', 'clingen_null'] → column: vg_predicted_perm
```

Join on `(source_mapped, gene_symbol)` where `background_null → background`, `clingen_null → clingen`.

### Per-region Vg depletion (obs / perm)

```
obs  = all_genes where source in ['background', 'clingen']          → columns: vg_{region}
perm = all_genes where source in ['background_null', 'clingen_null'] → columns: vg_{region}_perm
```

Same join on `(source_mapped, gene_symbol)`.

---

## Wrong Patterns — Do NOT Do These

```python
# WRONG: using _perm from the real gnomAD files as the null expectation
ratio = df_real['vg_predicted'] / df_real['vg_predicted_perm']

# WRONG: using vg_{region}_perm from the real files for per-region depletion
ratio = df_real['vg_promoter_core'] / df_real['vg_promoter_core_perm']

# WRONG: using vg_{region} (non-perm) from null files — these are zero
perm = df_null['vg_promoter_core']  # always 0, useless
```

---

## Rationale

The null files (`_null`) use a **different gene-window design** — AFs are permuted across the genome while keeping the same variant set. The `vg_{region}_perm` columns in those files hold Vg recalculated under permuted AFs and represent the sequence-context null expectation. The `_perm` columns inside the real files (`background`, `clingen`) exist as pipeline byproducts and do not represent the intended null model for purifying selection analyses.
