# Getting libraries
library("limma")
library("ggplot2")
library("tidyverse")
library("MASS")

# Setting working directory (line will be missing in script)
setwd("/projects/kaczorowski-lab/USERS/saulm/adbxd_proteomics/")

# Getting gemma_bxd() function
source("./code/gemma_bxd.R")

# Loading data
adbxd_hippo_pgroup            = readRDS("./data/adbxd_hippo_pgroup_annotated.RDS")
adbxd_hippo_sampanno          = readRDS("./data/adbxd_hippo_sampanno_quantnorm.RDS")
adbxd_genotypes               = readRDS("./data/bxd_geno_data.RDS")
adbxd_geno_meta               = readRDS("./data/bxd_geno_meta.RDS")

# Getting rank-normal transform function
norm_rank_transform = function(x, c = 0) {
  stopifnot(is.numeric(x) & is.vector(x))
  x_noNA = which(!is.na(x))
  N = length(x_noNA)
  x[x_noNA] = qnorm((rank(x[x_noNA], ties.method = "average") - c) / (N - (2 * c) + 1))
  return(x)
}

# Setting job array information
# Will adaptively change depending upon how many jobs are in array
# so long as array is formatted 1-X (starts at 1 and goes to a large number)
slurm_array_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
slurm_njobs = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX"))

# Getting information needed to subset the data.

# Getting rows vector
nrow_pgroup = nrow(adbxd_hippo_pgroup)
nrows_in_vector = ceiling(nrow_pgroup / slurm_njobs)
rows_vector = (1 + (slurm_array_id * nrows_in_vector - nrows_in_vector)):ifelse(slurm_array_id != slurm_njobs,
                                                                                (slurm_array_id * nrows_in_vector),
                                                                                nrow(adbxd_hippo_pgroup))

# Making genotype matrix
adbxd_genodata = cbind(data.frame(id = row.names(adbxd_genotypes),
                                  x = rep("X", times = nrow(adbxd_genotypes)),
                                  y = rep("Y", times = nrow(adbxd_genotypes))),
                       adbxd_genotypes)

# Making annotation
adbxd_geno_annot = adbxd_geno_meta
adbxd_geno_annot["14P_no_data","Mb_mm39"] = 3.000000
adbxd_geno_annot$id = row.names(adbxd_geno_annot)
adbxd_geno_annot$pos = adbxd_geno_annot$Mb_mm39 * 1e6
adbxd_geno_annot = adbxd_geno_annot[,c("id","pos","Chr")]

# Setting minimum number of rows to be present for analysis
min_nrow = 10

for (i in rows_vector) {
  protein_i = row.names(adbxd_hippo_pgroup)[i]
  adbxd_pabund_i = adbxd_hippo_sampanno
  adbxd_pabund_i$protein_i = adbxd_hippo_pgroup[protein_i,adbxd_pabund_i$sample_id]
  adbxd_mapping_i = adbxd_pabund_i |>
    filter(sample_type == "Hippocampus individual" &
             !(strain %in% c("C57BL/6J","DBA/2J"))   &
             genotype == "Ntg" &
             !is.na(protein_i) &
             age == "Months06") |>
    group_by(strain) |>
    summarize(mean = mean(protein_i),
              sd = sd(protein_i),
              n = length(protein_i)) |>
    mutate(se = sd / sqrt(n))
  # Only mapping when at least min_nrow values are present.
  if (nrow(adbxd_mapping_i) >= min_nrow) {
    mapping_result_i = gemma_bxd(pheno = pull(adbxd_mapping_i,"mean"),
                                 geno  = adbxd_genodata[,c("id","x","y",pull(adbxd_mapping_i, "strain"))],
                                 anno  = adbxd_geno_annot,
                                 loco  = TRUE,
                                 perms = 0)
    mapping_result_i$protein = rep(protein_i, times = nrow(mapping_result_i))
    if (!exists("adbxd_pqtl_ntg_months06")) {
      adbxd_pqtl_ntg_months06 = mapping_result_i
    } else {
      adbxd_pqtl_ntg_months06 = rbind(adbxd_pqtl_ntg_months06,
                                      mapping_result_i)
    }
  }
}

for (i in rows_vector) {
  protein_i = row.names(adbxd_hippo_pgroup)[i]
  adbxd_pabund_i = adbxd_hippo_sampanno
  adbxd_pabund_i$protein_i = adbxd_hippo_pgroup[protein_i,adbxd_pabund_i$sample_id]
  adbxd_mapping_i = adbxd_pabund_i |>
    filter(sample_type == "Hippocampus individual" &
             !(strain %in% c("C57BL/6J","DBA/2J"))   &
             genotype == "5xFAD" &
             !is.na(protein_i) &
             age == "Months06") |>
    group_by(strain) |>
    summarize(mean = mean(protein_i),
              sd = sd(protein_i),
              n = length(protein_i)) |>
    mutate(se = sd / sqrt(n))
  # Only mapping when at least min_nrow values are present.
  if (nrow(adbxd_mapping_i) >= min_nrow) {
    mapping_result_i = gemma_bxd(pheno = pull(adbxd_mapping_i,"mean"),
                                 geno  = adbxd_genodata[,c("id","x","y",pull(adbxd_mapping_i, "strain"))],
                                 anno  = adbxd_geno_annot,
                                 loco  = TRUE,
                                 perms = 0)
    mapping_result_i$protein = rep(protein_i, times = nrow(mapping_result_i))
    if (!exists("adbxd_pqtl_5xfad_months06")) {
      adbxd_pqtl_5xfad_months06 = mapping_result_i
    } else {
      adbxd_pqtl_5xfad_months06 = rbind(adbxd_pqtl_5xfad_months06,
                                      mapping_result_i)
    }
  }
}

for (i in rows_vector) {
  protein_i = row.names(adbxd_hippo_pgroup)[i]
  adbxd_pabund_i = adbxd_hippo_sampanno
  adbxd_pabund_i$protein_i = adbxd_hippo_pgroup[protein_i,adbxd_pabund_i$sample_id]
  adbxd_mapping_i = adbxd_pabund_i |>
    filter(sample_type == "Hippocampus individual" &
             !(strain %in% c("C57BL/6J","DBA/2J"))   &
             genotype == "Ntg" &
             !is.na(protein_i) &
             age == "Months14") |>
    group_by(strain) |>
    summarize(mean = mean(protein_i),
              sd = sd(protein_i),
              n = length(protein_i)) |>
    mutate(se = sd / sqrt(n))
  # Only mapping when at least min_nrow values are present.
  if (nrow(adbxd_mapping_i) >= min_nrow) {
    mapping_result_i = gemma_bxd(pheno = pull(adbxd_mapping_i,"mean"),
                                 geno  = adbxd_genodata[,c("id","x","y",pull(adbxd_mapping_i, "strain"))],
                                 anno  = adbxd_geno_annot,
                                 loco  = TRUE,
                                 perms = 0)
    mapping_result_i$protein = rep(protein_i, times = nrow(mapping_result_i))
    if (!exists("adbxd_pqtl_ntg_months14")) {
      adbxd_pqtl_ntg_months14 = mapping_result_i
    } else {
      adbxd_pqtl_ntg_months14 = rbind(adbxd_pqtl_ntg_months14,
                                      mapping_result_i)
    }
  }
}

for (i in rows_vector) {
  protein_i = row.names(adbxd_hippo_pgroup)[i]
  adbxd_pabund_i = adbxd_hippo_sampanno
  adbxd_pabund_i$protein_i = adbxd_hippo_pgroup[protein_i,adbxd_pabund_i$sample_id]
  adbxd_mapping_i = adbxd_pabund_i |>
    filter(sample_type == "Hippocampus individual" &
             !(strain %in% c("C57BL/6J","DBA/2J"))   &
             genotype == "5xFAD" &
             !is.na(protein_i) &
             age == "Months14") |>
    group_by(strain) |>
    summarize(mean = mean(protein_i),
              sd = sd(protein_i),
              n = length(protein_i)) |>
    mutate(se = sd / sqrt(n))
  # Only mapping when at least min_nrow values are present.
  if (nrow(adbxd_mapping_i) >= min_nrow) {
    mapping_result_i = gemma_bxd(pheno = pull(adbxd_mapping_i,"mean"),
                                 geno  = adbxd_genodata[,c("id","x","y",pull(adbxd_mapping_i, "strain"))],
                                 anno  = adbxd_geno_annot,
                                 loco  = TRUE,
                                 perms = 0)
    mapping_result_i$protein = rep(protein_i, times = nrow(mapping_result_i))
    if (!exists("adbxd_pqtl_5xfad_months14")) {
      adbxd_pqtl_5xfad_months14 = mapping_result_i
    } else {
      adbxd_pqtl_5xfad_months14 = rbind(adbxd_pqtl_5xfad_months14,
                                        mapping_result_i)
    }
  }
}

saveRDS(adbxd_pqtl_ntg_months06,
        paste0("./data/pqtl/new_pqtl/adbxd_pqtl_ntg_months06_a",slurm_array_id,".RDS"))
saveRDS(adbxd_pqtl_5xfad_months06,
        paste0("./data/pqtl/new_pqtl/adbxd_pqtl_5xfad_months06_a",slurm_array_id,".RDS"))
saveRDS(adbxd_pqtl_ntg_months14,
        paste0("./data/pqtl/new_pqtl/adbxd_pqtl_ntg_months14_a",slurm_array_id,".RDS"))
saveRDS(adbxd_pqtl_5xfad_months14,
        paste0("./data/pqtl/new_pqtl/adbxd_pqtl_5xfad_months14_a",slurm_array_id,".RDS"))
