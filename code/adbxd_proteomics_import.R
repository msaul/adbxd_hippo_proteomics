library("readxl")
library("janitor")

adbxd_metadata = read_xlsx(path = "~/Downloads/adbxd_proteomics/MouAD-Contract2-2021-Hipp-MetaData-All.xlsx",
                           sheet = "Hipp-Meta-UW-MatchC1") |>
  clean_names()
adbxd_metadata = as.data.frame(adbxd_metadata)
row.names(adbxd_metadata) = gsub("-","_",adbxd_metadata$uw_sample_id)

adbxd_metadata_jax = read_xlsx(path = "~/Downloads/adbxd_proteomics/MouAD-Contract2-2021-Hipp-MetaData-All.xlsx",
                               sheet = "Hipp-Meta-Jax") |>
  clean_names()

adbxd_metadata_uw = read_xlsx(path = "~/Downloads/adbxd_proteomics/MouAD-Contract2-2021-Hipp-MetaData-All.xlsx",
                              sheet = "Hipp-Meta-UW") |>
  clean_names()

adbxd_proteome = read.csv("~/Downloads/adbxd_proteomics/20220504-MouAD-Hipp-C2-protGrp_batchadj.csv")
adbxd_annot = data.frame(row.names = adbxd_proteome$Protein,
                         protein_id = adbxd_proteome$Protein,
                         n_ids = lengths(strsplit(adbxd_proteome$Protein, ",")),
                         prefix = rep("", times = nrow(adbxd_proteome)),
                         uniprot = rep("", times = nrow(adbxd_proteome)),
                         id = rep("", times = nrow(adbxd_proteome)),
                         peptide = adbxd_proteome$Peptide)
adbxd_proteome = as.matrix(adbxd_proteome[,-1:-2])
colnames(adbxd_proteome) = gsub("\\.","_",colnames(adbxd_proteome))
row.names(adbxd_proteome) = row.names(adbxd_annot)

adbxd_metadata = adbxd_metadata[colnames(adbxd_proteome),]


for (i in 1:nrow(adbxd_annot)) {
  protein_ids_i = unlist(strsplit(adbxd_annot[i,"protein_id"],","))
  for (j in 1:adbxd_annot[i,"n_ids"]) {
    if (j == 1) {
      prefix_i = gsub("^\\s+(.*)\\|(.*)\\|(.*)$","\\1",protein_ids_i[j])
      uniprot_i = gsub("^\\s+(.*)\\|(.*)\\|(.*)$","\\2",protein_ids_i[j])
      id_i = gsub("^\\s+(.*)\\|(.*)\\|(.*)$","\\3",protein_ids_i[j])
    } else {
      prefix_i = c(prefix_i, gsub("^\\s+(.*)\\|(.*)\\|(.*)$","\\1",protein_ids_i[j]))
      uniprot_i = c(uniprot_i, gsub("^\\s+(.*)\\|(.*)\\|(.*)$","\\2",protein_ids_i[j]))
      id_i = c(id_i, gsub("^\\s+(.*)\\|(.*)\\|(.*)$","\\3",protein_ids_i[j]))
    }
  }
  adbxd_annot[i,"prefix"] = paste(unique(prefix_i), sep = "", collapse = ";")
  adbxd_annot[i,"uniprot"] = paste(unique(uniprot_i), sep = "", collapse = ";")
  adbxd_annot[i,"id"] = paste(unique(id_i), sep = "", collapse = ";")
}

# Getting all Swiss-Prot annotated data.
adbxd_annot_filtered = adbxd_annot[which(adbxd_annot$n_ids == 1 & adbxd_annot$prefix == "sp"),]
n_queries = 1 + nrow(adbxd_annot_filtered) %/% 100

for (i in 1:n_queries) {
  if (i != n_queries) {
    indices_i = ((i - 1) * 100 + 1):(i * 100)
  } else {
    indices_i = ((i - 1) * 100 + 1):(nrow(adbxd_annot_filtered))
  }
  query_i = adbxd_annot_filtered$uniprot[indices_i]
  requestURL_i = paste0("https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession=",
                        paste(query_i, sep = "", collapse = "%2C"))
  r_i = httr::GET(requestURL_i, httr::accept("application/json"))

  httr::stop_for_status(r_i)

  json_i = httr::content(r_i)
  # names(json_i) = query_i
  # df_i = jsonlite::fromJSON(json_i, flatten = FALSE)
  # row.names(df_i) = query_i

  if (i == 1) {
    annot_list = json_i
  } else {
    annot_list = c(annot_list, json_i)
  }

  cat("Iteration ", i, "\n", sep = "")
}

for (i in 1:n_queries) {
  if (i != n_queries) {
    indices_i = ((i - 1) * 100 + 1):(i * 100)
  } else {
    indices_i = ((i - 1) * 100 + 1):(nrow(adbxd_annot_filtered))
  }
  query_i = adbxd_annot_filtered$uniprot[indices_i]
  requestURL_i = paste0("https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession=",
                        paste(query_i, sep = "", collapse = "%2C"))
  r_i = httr::GET(requestURL_i, httr::accept("application/json"))

  httr::stop_for_status(r_i)

  json_i = httr::content(r_i)
  # names(json_i) = query_i
  # df_i = jsonlite::fromJSON(json_i, flatten = FALSE)
  # row.names(df_i) = query_i

  if (i == 1) {
    annot_list_2 = json_i
  } else {
    annot_list_2 = c(annot_list_2, json_i)
  }

  cat("Iteration ", i, "\n", sep = "")
}

names(annot_list) = unlist(lapply(annot_list, function(X) {return(X[["accession"]])}))
names(annot_list_2) = unlist(lapply(annot_list_2, function(X) {return(X[["accession"]])}))

chr_list = unlist(lapply(annot_list, function(X){return(X[["gnCoordinate"]][[1]]$genomicLocation$chromosome)}))

