#' Peptide/Protein group inference
#'
#' @param dat input dataframe
#' @param protGrp logical variable indicating if protein groups should be inferred
#' @param average logical variable indicating if abundances should be averaged
#'
#' @return dataframe with the peptide/protein groups inferred and abundances averaged or summed
#' @import stats
#' @import utils
#' @importFrom data.table .SD
#' @export
#'
#' @examples
#' pepnames <- c("a@A", "b@A", "b@B", "c@C", "d@C", "a@D", "b@D")
#' dat <- matrix(runif(35, min=8, max=25), 7, 5)
#' rownames(dat) <- pepnames
#' pepGrp <- get_pgroup(dat)

get_pgroup <- function(dat, protGrp = F, average = T) {

  # warnings off
  w <- getOption("warn")
  options(warn = -1)

  # global variable
  ..i <- NULL

  # extract peptide id
  pep = sapply(strsplit(row.names(dat),"@"), "[", 1)
  # extract protein id
  prot = sapply(strsplit(row.names(dat),"@"), "[", 2)
  prot = sapply(strsplit(prot,"_"), "[", 1)
  # create dataframe with peptide and protein ids populated
  pep_prot = cbind(pep,prot,dat)

  # group proteins
  if(protGrp == T) {
    pep_prot = data.table::setDT(as.data.frame(
      pep_prot))[, lapply(.SD, function(x) toString(na.omit(x))), by = prot]
  }

  # group peptides
  pep_prot = data.table::setDT(as.data.frame(
    pep_prot))[, lapply(.SD, function(x) toString(na.omit(x))), by = pep]

  # initialize group sum and average abundances
  prot.grp.avg = pep_prot[,1:2]
  prot.grp.sum = pep_prot[,1:2]

  # protein group sum and average abundances
  if (protGrp == T) {
    if (average == T) {
      # calculate protein group average
      for (i in 3:(ncol(dat)+2)) {
          avg = sapply(as.matrix(pep_prot[, ..i]),
                       function(x) mean(scan(text = x, what=numeric(),sep=",",
                                           quiet = T), na.rm=TRUE))
        prot.grp.avg = cbind(prot.grp.avg, avg)
      }
      colnames(prot.grp.avg) = colnames(pep_prot)
    } else {
      # calculate protein group sum
      for (i in 3:(ncol(dat)+2)) {
          sum = sapply(as.matrix(pep_prot[, ..i]),
                       function(x) sum(scan(text = x, what=numeric(),sep=",",
                                            quiet = T), na.rm=TRUE))
        prot.grp.sum = cbind(prot.grp.sum, sum)
      }
      colnames(prot.grp.sum) = colnames(pep_prot)
    }
  } else {
    if (average == T){
      # subset the first abundance (same as average)
      for (i in 3:(ncol(dat)+2)){
      avg = sapply(strsplit(pep_prot[[i]], ","), "[", 1)
      prot.grp.avg = cbind(prot.grp.avg, avg)
    }
    colnames(prot.grp.avg) = colnames(pep_prot)
    } else {
      # sum peptide group abundances
        for (i in 3:(ncol(dat)+2)) {
            sum = sapply(as.matrix(pep_prot[, ..i]),
                         function(x) sum(scan(text = x, what=numeric(),sep=",",
                                              quiet = T), na.rm=TRUE))
          prot.grp.sum = cbind(prot.grp.sum, sum)
        }
        colnames(prot.grp.sum) = colnames(pep_prot)
      }
      prot.grp = prot.grp.sum
    }

  # return averaged or summed abundances based on average parameter
  if (average == T)
    prot.grp = prot.grp.avg
  else
    prot.grp = prot.grp.sum

  # drop the first two columns of the dataframe containing ids
  prot.grp = as.matrix(prot.grp[,-c(1,2)])
  # change class of the dataframe to numeric column-wise
  prot.grp = apply(prot.grp, 2, as.numeric)
  # set row names with new ids
  rownames(prot.grp) = paste0(pep_prot$pep, " @ ", pep_prot$prot)

  # warnings back on
  options(warn = w)

  # return peptide/protein grouped abundance data matrix
  return(prot.grp)
}

