## Directory Structure

---

- /data/
  - chromatin_associated_genes/
    - &lt;gene_name&gt;/
      - ss_final_out_&lt;GENE&gt;/
      - &lt;GENE&gt;_cv_glmnet_metrics_tumor.txt
      - &lt;GENE&gt;_enrich_res/
  - survival/
    - curated_survival_BRCA.txt
    - brca_metadata/
  - tissues/
    - tcga_metadata/
      - &lt;tissue&gt;/
        - bam/
  - counts/
    - exonLevel_allCounts_filtbyexpr.tsv
  - salmon/
    - salmon_predictors/
      - salmon_pvt1tpm_fix.tsv

---

## Noticeable File Path Examples

- **Splicing efficiency dataframe (PVT1)**:  
  `/data/chromatin_associated_genes/pvt1/condition_altSS_concat_Control_Tumor_SE_SS3.altSS_L50.tsv`  
  _(similar for other genes, personalized for each)_

- **CV Elastic Net metrics**:  
  `/data/chromatin_associated_genes/<gene>/<gene>_cv_glmnet_metrics_tumor.txt`

- **Survival data**:  
  `/data/survival/brca_metadata/curated_survival_BRCA.txt`

- **Exon counts**:  
  `/data/counts/exonLevel_allCounts_filtbyexpr.tsv`  
  _or `allExon_allSS_RPM.txt` for 670×30887 genes + 34 ss (includes extra columns)_

- **Salmon 19 TPM dataframe**:  
  `/data/salmon/salmon_predictors/salmon_pvt1tpm_Subfix.tsv`

- **TCGA manifest metadata**:  
  `/data/tissues/tcga_metadata/<tissue>/bam/`  
  _BRCA manifest: `/data/concat_manifest_with_details_V4.tsv` (with full 2382 miRNA targets)_

- **Generalization results**:  
  `/data/tissues/<tissue>/PVT1_generalized_glmnet_coeff_metrics.tsv`



k-means code : https://github.com/evgnt/chrTT-seq/
Random Forest code: 
