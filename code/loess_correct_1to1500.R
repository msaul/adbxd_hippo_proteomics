# Loess correct: 1 to 1500
# Michael C. Saul
# michael.saul [at] jax.org
# 2022-10-06

# This script batch corrects the proteomics data in parallel.
# It splits the proteomics data matrix into 1,500 smaller chunks.
# These chunks run in parallel on the cluster.

# The script makes use of a column index and a means of producing approximately
# equal row indices.

# Getting working directory
setwd("/projects/kaczorowski-lab/USERS/saulm/adbxd_proteomics")

# Loading libraries
library("tidyverse")
library("tibble")
library("proBatch")

# Loading the data
adbxd_hippo_quant = readRDS("./data/adbxd_hippo_peptide_quantnorm.RDS")
adbxd_hippo_sampanno = readRDS("./data/adbxd_hippo_sampanno_quantnorm.RDS")

# Getting SLURM ID
slurm_array_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Getting information needed to subset the data.

# Getting rows vector
nrow_quant = nrow(adbxd_hippo_quant)
rows_slurm_id = ((slurm_array_id - 1) %/% 15) + 1
nrows_in_vector = ceiling(nrow_quant / 100)
rows_vector = (1 + (rows_slurm_id * nrows_in_vector - nrows_in_vector)):ifelse(rows_slurm_id != 100, (rows_slurm_id * nrows_in_vector), nrow_quant)

# Getting col vector
cols_slurm_id = ((slurm_array_id - 1)  %% 15) + 1

# Getting matrix slice
adbxd_hippo_sampanno_subset = adbxd_hippo_sampanno[which(adbxd_hippo_sampanno$loess_correction_n == cols_slurm_id),]
adbxd_hippo_sampanno_subset_tibble = as_tibble(adbxd_hippo_sampanno_subset)
adbxd_hippo_quant_long = matrix_to_long(adbxd_hippo_quant[rows_vector,adbxd_hippo_sampanno$sample_id])

# Running Loess regression
adbxd_hippo_loess_long = adjust_batch_trend_df(adbxd_hippo_quant_long,
                                               sample_annotation = adbxd_hippo_sampanno_subset_tibble,
                                               batch_col = "batch",
                                               order_col = "order")

adbxd_hippo_loess = long_to_matrix(adbxd_hippo_loess_long)

# Saving out data
saveRDS(adbxd_hippo_loess_long,paste0("./data/loess/adbxd_hippo_loess_long_array_",slurm_array_id,".RDS"))
saveRDS(adbxd_hippo_loess,paste0("./data/loess/adbxd_hippo_loess_array_",slurm_array_id,".RDS"))
