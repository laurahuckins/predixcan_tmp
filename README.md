# predixcan_tmp

brief overview of files to create predixcan models.

Requirements:
1. RNAseq data (raw counts)
2. Covariates (PCA/Ancestry EVs)
3. Dosage data

Steps to running glmnet:

1. QC: SNPS, GENES, INDIVS (see QC files in this folder)
2. Run elastic net regresssion for each gene (see run_glmnet.R)
3. QC results for (2): see check_results.sh
4. Make Database (See xxx.py)
5. Make covariance matrices (See xxxxx)
