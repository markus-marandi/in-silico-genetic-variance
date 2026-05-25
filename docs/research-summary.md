# Research Summary

This repository implements the computational pipeline behind the thesis **Deep-Learning Prediction of Human Gene Expression Variability from Genetic Variation**. It is designed to turn population and synthetic variant sets into gene-level estimates of predicted cis-regulatory genetic variance (`V_G`) for downstream statistical analysis.

## Thesis Context

The thesis studies dosage-sensitive genes, with a focus on ClinGen haploinsufficient genes. The motivating gap is that established constraint metrics such as LOEUF, pLI, pTriplo, ncGERP, ncCADD, and Episcore quantify coding intolerance, structural variation, conservation, or epigenomic annotations, but do not directly measure predicted depletion of cis-regulatory expression variance.

The pipeline adapts the ANEVA-style genetic-variance formulation to deep-learning variant-effect predictions:

```text
V_G = sum_i 2 p_i (1 - p_i) beta_i^2
```

Here, `p_i` is the allele frequency of variant `i`, and `beta_i` is the AlphaGenome-predicted expression effect size for that variant.

## Study Design Implemented by the Pipeline

The thesis analysis used two gene sets:

- 316 ClinGen haploinsufficient genes classified as having sufficient evidence for haploinsufficiency.
- 349 background coding genes sampled from protein-coding genes in the gnomAD v4.1 constraint table, excluding high-confidence ClinGen haploinsufficient genes.

Both sets were filtered for MANE Select transcript availability, skeletal-muscle expression, and promoter non-overlap. For each retained gene, the pipeline considers variants within +/-10 kb of the transcription start site on GRCh38.

## Variant Sets

### Observed gnomAD Variants

Observed variants are extracted from gnomAD v4.1 genomes within each TSS-centered window. No allele-frequency threshold is applied, so the observed set spans singletons, rare variants, and common variants.

### Synthetic Null Variants

The synthetic null starts from all possible SNVs in the same gene windows. Variants observed in gnomAD are removed by anti-join. The remaining possible SNVs are sampled to match the observed mutation spectrum and CpG methylation context, then downsampled per gene to the observed variant count.

The synthetic null provides a counterfactual mutational baseline, but it does not fully isolate purifying selection from sequence composition. In the thesis, ClinGen haploinsufficient promoters remained lower in null `V_G`, consistent with CpG-rich promoter architecture influencing both mutational opportunity and predicted effect sizes.

## AlphaGenome Scoring

AlphaGenome is used as a sequence-to-expression model. For each variant, the model compares reference and alternate alleles in a 1 Mb sequence context, while the thesis variant set itself is restricted to +/-10 kb around the TSS. RNA-seq signal is aggregated over annotated exonic bins of the target gene, yielding a gene-level log2 fold-change effect size.

The primary thesis signal is skeletal-muscle polyA+ RNA-seq, selected to match the GTEx skeletal-muscle expression and eQTL validation context.

## Aggregation Contract

The prepare-results stage keeps RNA protocols explicit:

- `track_key` identifies an exact AlphaGenome RNA track.
- `protocol_group` groups related RNA protocols such as `polyA_plus_rna_seq` and `total_rna_seq`.

Variant-level outputs may therefore contain multiple rows for the same `gene_id` and `variant_id` when distinct RNA tracks are present. Deduplication and aggregation must be performed within the appropriate track/protocol identity, not by collapsing all RNA rows together.

For thesis observed-versus-null analyses:

| Dataset type | Use for total `V_G` | Use for AF partitions | Use for regional windows |
| --- | --- | --- | --- |
| Observed gnomAD | `vg_predicted` | `vg_common`, `vg_rare` | `vg_*` |
| Synthetic/null | `vg_predicted_perm` | `vg_common_perm`, `vg_rare_perm` | `vg_*_perm` |

Permutation or sanity columns in real files are quality-control artifacts. They should not be interpreted as the biological null model.

## Interpretation

The downstream thesis analysis found that gene-level predicted `V_G` correlated with empirical eQTL-derived `V_G`, whereas matched per-variant eQTL concordance was not supported. ClinGen haploinsufficient genes carried lower total predicted regulatory variance than background genes, but null-model and spatial-architecture results require careful interpretation because sequence composition remains a confounder.

Predicted `V_G` should be treated as a gene-level research and prioritization metric, not as clinical evidence for individual variants.

## Acknowledgements and Funding

This work was supervised by Philipp Rentzsch and Tuuli Lappalainen at the Lappalainen Lab, KTH / Science for Life Laboratory.

Computational resources were provided by NAISS at PDC Center for High Performance Computing, KTH Royal Institute of Technology, allocation `NAISS 2024/6-322`, and at NSC, Linköping University, allocation `Berzelius-2025-176`.

Funding support came from a Wallenberg Scholar award 2024 to T. Lappalainen, Knut and Alice Wallenberg Foundation grant `KAW 2023.0337`, and a Göran Gustavsson award 2023 to T. Lappalainen.
