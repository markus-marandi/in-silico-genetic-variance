# In-Silico Genetic Variance Pipeline

This repository contains a computational pipeline designed to predict gene expression using AlphaGenome. The workflow generates variant tables from gene lists using gnomad data and additionally also null variant sets. Then they are scored  for molecular impact, and aggregated to variant-level predictions and into gene-level summaries. 

## Pipeline Overview

The workflow is divided into three distinct stages:
1.  **Generation:** Creating variant sets (gnomAD/Null).
2.  **Scoring:** Scoring via the AlphaGenome API.
3.  **Aggregation:** Calculating $V_{G,pred}$ and gene-level metrics.

---

### Stage 1: Variant Generation

This stage prepares the input data for the model. It handles the generation of background variants, sampling from gnomAD, and formatting genomic intervals.

![Workflow Diagram](docs/figs/Figure%201.9.1:%20Computational%20workflow%20for%20variant%20generation%20and%20scoring.png)

**Key Objectives:**
* Select target genes (MANE Select transcripts).
* Generate SNPs within a defined window (1Mb around TSS).
* Annotate variants with observed Allele Frequencies (AF) where available.

---

### Stage 2: AlphaGenome Scoring

This stage submits the formatted variants to the AlphaGenome API to predict their molecular impact (log-fold change). This step is time-consuming, but smaller batch sizes do not require large amount of computational power. 

**Configuration Required (.env)**
To run this stage, you must authenticate with the AlphaGenome API. You need to create a `.env` file in the `2_score_variants_with_alphagenome/` directory containing your API key.

**Setup:**
1.  Create the file: `2_score_variants_with_alphagenome/.env`
2.  Add your key:
    ```bash
    API_KEY_PERSONAL=your_alphagenome_api_key_here
    ```

**Execution:**
The scoring logic handles token management, batching, and retries automatically.
```bash
# Example SLURM submission
sbatch --export=DATASET_ID=dataset4,SAMPLE_ID=background ... run_scoring.sh

```

---

### Stage 3: Aggregation & Analysis

The final stage processes the raw prediction chunks. It stitches them into a unified dataset, backfills metadata, and performs the mathematical reduction to gene-level metrics.

**Predicted Genetic Variance**
For every gene, we calculate the cumulative variance explained by the model using the formula:

$$ V_{G,pred} = \sum_{i \in \text{variants}} 2 \cdot p_i \cdot (1 - p_i) \cdot \beta_i^2 $$

Where:

* : Allele Frequency (AF)
* : Predicted effect size (`raw_score`)

**Outputs:**

* **Variant-Level Parquet:** Annotated list of all scored variants.
* **Gene-Level Parquet:** Aggregated metrics including , mean absolute effects, and spatial effect distributions (Promoter vs. Gene Body).

## Installation & Environment

This pipeline relies on `polars` for high-performance data manipulation.

```bash
# Install dependencies
pip install -r requirements.txt

```

Ensure your `sys.path` is correctly configured if running scripts directly from the `modules/` folder, or use the provided runner scripts which handle imports automatically.