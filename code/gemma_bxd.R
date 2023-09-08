# This is a wrapper function for GEMMA
# It is designed to simplify calling GEMMA from R and returning data from GEMMA

gemma_bxd = function(pheno,
                     geno,
                     anno,
                     covar = NULL,
                     gxe = NULL,
                     loco = TRUE,
                     save_raw = FALSE,
                     perms = 0) {
  # Throwing errors for some specific anticipated issues
  if (loco & !is.null(gxe)) {
    stop("GxE is currently incompatible with LOCO")
  } else if (is.null(anno) & !is.null(loco)) {
    stop("An annotation must be provided for LOCO")
  } else if (!is.data.frame(geno)) {
    stop("Genotype data are not formatted as a data.frame")
  } else if (length(perms) != 1 | !is.numeric(perms) | perms < 0) {
    stop("perms argument is not a single number greater than or equal to 0")
  }

  # Making temporary directory for data
  gemma_dir = tempdir(check = TRUE)
  gemma_outdir = paste0(gemma_dir,"/output"); dir.create(gemma_outdir)

  # Writing out annotation
  write.table(anno,
              paste0(gemma_dir,"/adata.txt"),
              sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)

  # Writing out genotype and phenotype data
  write.table(pheno,
              paste0(gemma_dir,"/pdata.txt"),
              sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  write.table(geno,
              paste0(gemma_dir,"/gdata.txt"),
              sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)

  # Making commands
  gemma_kin_command = paste0("gemma -g ",
                             gemma_dir, "/gdata.txt -p ",
                             gemma_dir, "/pdata.txt",
                             ifelse(!is.null(anno), paste0(" -a ", gemma_dir, "/adata.txt"), ""))
  gemma_gwa_command = paste0("gemma -g ",
                             gemma_dir, "/gdata.txt -p ",
                             gemma_dir, "/pdata.txt",
                             ifelse(!is.null(anno), paste0(" -a ", gemma_dir, "/adata.txt"), ""))

  # Making covariance data
  if (!is.null(covar)) {
    write.table(covar,
                paste0(gemma_dir,"/cdata.txt"),
                sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    gemma_kin_command = paste0(gemma_kin_command, " -c ", gemma_dir, "/cdata.txt")
    gemma_gwa_command = paste0(gemma_gwa_command, " -c ", gemma_dir, "/cdata.txt")
  }

  # Making GxE data
  if (!is.null(gxe)) {
    # stop("GxE doesn't currently work") # GxE isn't working for now.
    write.table(gxe,
                paste0(gemma_dir,"/gxedata.txt"),
                sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    gemma_kin_command = paste0(gemma_kin_command, " -gxe ", gemma_dir, "/gxedata.txt")
    gemma_gwa_command = paste0(gemma_gwa_command, " -gxe ", gemma_dir, "/gxedata.txt")
  }

  # Making permutations files
  if (perms > 0) {
    for (i in 1:perms) {
      write.table(sample(pheno),
                  file = paste0(gemma_dir, "/pdata_perm",i,".txt"),
                  sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }

  # Running LOCO
  if (loco) {
    for (i in c(1:19,"X")) {
      system(paste0("cd ", gemma_dir, " ; ", gemma_kin_command, " -gk 1 -loco ",
                    i," -o kinship_chr",i,"_loco"),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
      # cat("Ran kinship for chr",i,"\n",sep="")
      system(paste0("cd ", gemma_dir, " ; ", gemma_gwa_command, " -k ",
                    gemma_outdir, "/kinship_chr",
                    i,"_loco.cXX.txt -lmm 2 -loco ", i,
                    " -lmin 0.01 -lmax 100 -o gwa_chr",i,"_loco"),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
      # cat("Ran GEMMA for chr",i,"\n",sep="")
      file_i = read.table(paste0(gemma_outdir,"/gwa_chr",i,"_loco.assoc.txt"),
                          sep = "\t", header = TRUE)
      file_i$lod = -1 * log10(file_i$p_lrt)
      if (i == "1") {
        gemma_results = file_i
      } else {
        gemma_results = rbind(gemma_results, file_i)
      }
      if (perms == 0 | !(exists("perms"))) {
        next
      } else {
        for (j in 1:perms) {
          permu_gemma_command = gsub("pdata\\.txt", paste0("pdata_perm",j,".txt"), gemma_gwa_command)
          system(paste0("cd ", gemma_dir, " ; ", permu_gemma_command, " -k ",
                        gemma_outdir, "/kinship_chr",
                        i,"_loco.cXX.txt -lmm 2 -loco ", i,
                        " -lmin 0.01 -lmax 100 -o /perm",j,"_chr",i,"_loco"),
                 ignore.stdout = TRUE, ignore.stderr = TRUE)
          permu_ij = read.table(paste0(gemma_outdir,"/perm",j,"_chr",i,"_loco.assoc.txt"),
                                sep = "\t", header = TRUE)
          permu_ij$lod = -1 * log10(permu_ij$p_lrt)
          if (i == "1" & j == 1) {
            permus = data.frame(chr = i,
                                permu = j,
                                lod = max(permu_ij$lod))
          } else {
            permus = rbind(permus,
                           data.frame(chr = i,
                                      permu = j,
                                      lod = max(permu_ij$lod)))
          }
        }
      }
    }
  } else {
    system(paste0("cd ", gemma_dir, " ; ", gemma_kin_command, " -gk 1 -o kinship"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    # cat("Ran kinship\n", sep = "")
    system(paste0("cd ", gemma_dir, " ; ", gemma_gwa_command, " -k ", gemma_outdir,
                  "/kinship.cXX.txt -lmm 2 -lmin 0.01 -lmax 100 -o gwa"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    # cat("Ran GEMMA\n", sep = "")
    gemma_results = read.table(paste0(gemma_outdir,"/gwa.assoc.txt"),
                               sep = "\t", header = TRUE)
    gemma_results$lod = -1 * log10(gemma_results$p_lrt)
    if (exists("perms") & perms > 0) {
      for (j in 1:perms) {
        permu_gemma_command = gsub("pdata\\.txt", paste0("pdata_perm",j,".txt"), gemma_gwa_command)
        system(paste0("cd ", gemma_dir, " ; ", permu_gemma_command, " -k ", gemma_outdir,
                      "/kinship.cXX.txt -lmm 2 -lmin 0.01 -lmax 100 -o perm",j),
               ignore.stdout = TRUE, ignore.stderr = TRUE)
        permu_j = read.table(paste0(gemma_outdir,"/perm",j,".assoc.txt"),
                             sep = "\t", header = TRUE)
        permu_j = aggregate(p_lrt ~ chr, data = permu_j, FUN = function(x) {max(-1 * log10(x))})
        permu_j$permu = rep(j, times = nrow(permu_j))
        colnames(permu_j) = c("chr","lod","permu")
        permu_j = permu_j[,c("chr","permu","lod")]
        if (j == 1) {
          permus = permu_j
        } else {
          permus = rbind(permus, permu_j)
        }
      }
    }
  }
  if (save_raw) {
    # system(paste0("cd ", gemma_outdir, " ; tar cvfz gemma_results.tar.gz *.txt"))
  }
  if (perms > 0) {
    gemma_results = list(gemma_results = gemma_results[which(gemma_results$chr %in% c(1:19,"X")),],
                         gemma_perms = permus[which(permus$chr %in% c(1:19,"X")),])
  } else {
    gemma_results = gemma_results[which(gemma_results$chr %in% c(1:19,"X")),]
  }

  # Deleting data
  system(paste0("rm -r ", gemma_outdir),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("rm -r ", gemma_dir, "/*data*.txt"),
         ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Returning results
  return(gemma_results)
}

plot_gemma_bxd = function(gemma_results,
                          threshold = 0.63,
                          gap = 5e6,
                          add_permu = TRUE,
                          genome_build = NULL,
                          colors = c("#008FC0","#05396B")) {
  stopifnot(require("ggplot2"))
  if (add_permu) {
    maxpermu = aggregate(lod ~ permu, data = gemma_results[[2]], FUN = function(x) {max(x)})
    permu_threshold = quantile(maxpermu$lod, probs = threshold)[1]
    gemma_df = gemma_results[[1]]
  } else {
    gemma_df = gemma_results
  }
  gemma_df = gemma_df[which(gemma_df$chr %in% c(1:19,"X")),]

  midpoints = c()
  start = 0
  gemma_df$calcpos = rep(-1, times = nrow(gemma_df))
  gemma_df$stagger = rep("", times = nrow(gemma_df))
  for (i in c(1:19,"X")) {
    stagger = ifelse((which(c(1:19,"X") %in% i) %% 2) == 1, "one", "two")
    gemma_df[which(gemma_df$chr == i),"stagger"] = stagger
    maxpos = max(gemma_df[which(gemma_df$chr == i),"ps"])
    minpos = min(gemma_df[which(gemma_df$chr == i),"ps"])
    midpoints[i] = (0.5 * (maxpos + minpos)) + start
    gemma_df[which(gemma_df$chr == i),"calcpos"] = start + gemma_df[which(gemma_df$chr == i),"ps"]
    start = start + maxpos + gap
  }
  plot = ggplot(data = gemma_df,
                aes(x = calcpos,
                    y = lod,
                    color = stagger,
                    group = chr)) +
    geom_line(linewidth = 0.5) +
    theme_bw() +
    scale_x_continuous(breaks = midpoints,
                       labels = names(midpoints)) +
    scale_color_manual(values = colors) +
    theme(legend.position = "none",
          panel.grid = element_line(color = "#FFFFFF")) +
    ylab("LOD")
  if (add_permu) {
    plot = plot +
      geom_hline(yintercept = permu_threshold,
                 linetype = "dashed", color = "#333F48")
  }
  if (is.null(genome_build)) {
    plot = plot + xlab("Position")
  } else {
    plot = plot + xlab(paste0("Position (", genome_build, ")"))
  }
  return(plot)
}
