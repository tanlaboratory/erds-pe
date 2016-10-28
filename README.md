# erds-pe
Introduction

ERDS-pe is a tool designed to detect detect CNVs from whole-exome sequencing (WES) data. ERDS-pe employs principal component analysis to normalize WES data and incorporates RD and single-nucleotide variation information together as a hybrid signal into a paired hidden Markov model to infer CNVs from WES data. 


Installation

Pipeline commands
1. Calculate RPKM values for all samples.
python erds_pe rpkm --target TARGET --input INPUT --output OUTPUT

2. Merge RPKM values of all samples into a matrix.
python erds_pe.py merge_rpkm --rpkm_dir RPKM_DIR --target TARGET  [--output OUTPUT]

3. Normalization
  python erds_pe.py svd \
  --rpkm_matrix RPKM_matrix_toy_data.raw

4. CNV calling
  python erds_pe.py discover \
  --params params.txt \
  --datafile RPKM_matrix_toy_data.raw.SVD \
  --sample NA12878 \
  --vcf NA12878_exome.vcf.gz \
  --hetsnp True \
  --tagsnp Ture \
  --tagsnp_file tagSNP_hg19.txt \
  --output NA12878.svd.Het.Tag.cnv


Contact
renjie.tan@outlook.com
