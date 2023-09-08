library("rawrr")

raw_files_dir = "/fastscratch/saulm/raw/"
raw_files = list.files(path = raw_files_dir,
                       pattern = "\\.raw$")

for (i in 1:length(raw_files)) {
  file_i = raw_files[i]
  header_i = unlist(rawrr::readFileHeader(paste0(raw_files_dir, file_i)))
  df_header_i = as.data.frame(matrix(c(file_i, header_i), nrow = 1))
  colnames(df_header_i) = c("raw_file", names(header_i))
  if (i == 1) {
    raw_files_header = df_header_i
  } else {
    raw_files_header = rbind(raw_files_header,
                             df_header_i)
  }
}

write.table(raw_files_header, "/projects/kaczorowski-lab/USERS/saulm/adbxd_proteomics/data/raw/raw_chrlib_file_headers.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
