# Use SVA to account for cell composition in WGS data
# Based on tutorial at https://akhilesh362.wordpress.com/
# and the SVA manual https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

# Install and load packages
setUpWorkSpace = function() {
  source("https://bioconductor.org/biocLite.R")
  
  installIfNeededFromBioconductor = function(pkg){
    if(pkg %in% rownames(installed.packages()) == FALSE) {
      biocLite(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
  }
  
  installIfNeeded = function(pkg){
    if(pkg %in% rownames(installed.packages()) == FALSE) {
      install.packages(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
  }
  
  
  biopackages = c("sva", "limma", "BiSeq", "qvalue")
  lapply(biopackages, installIfNeededFromBioconductor)
  
  packages = c("MASS", "ggplot2", "isva")
  lapply(packages, installIfNeeded)
  
  if(!file.exists("readMethylationSamples.R")){
    cat("Error: source script readMethylationSamples.R was not found")
    cat("Are you in the correct working directory?")
    quit(save="no")
  }
  source("readMethylationSamples.R") # simplify reading of bismark files
}

setUpWorkSpace()

data.env = new.env() # data for the ongoing analysis

# Variables to change before running on cluster
root.dir =  "Y:/"
# root.dir =  "/mnt/research2/"

sample.file = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra.txt")
cov.folder  = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/cov_files")

default.min.beta.diff   = 0.01 # minimum required difference in methylation fraction
default.min.reads       = 5   # minimum number of reads in every sample
default.high.exp.filter = 0.1  # fraction of highest coverage loci to exclude from each sample


#' Load and filter data
#' 
#' Create the table for SVA input.  The data should be a matrix with features
#' in the rows and samples in the columns.
#' 
#' Filter the incoming data by readcount. Discard the CpGs with the top 10% of total read coverage in each sample; this should
#' help remove artefacts due to PCR duplication
#' 
#' @param  min.reads Only include CpGs for which there are this many reads in every sample
#' @return a matrix with CpGs in rows and samples in columns.
loadData = function(min.reads){
  cat("Excluding loci with top 10% of read coverage in each sample\n")
  cat("Excluding loci with less than", min.reads, "reads in any sample\n")
  
  APPLY_TO_ROWS = 1 # parameter for base::apply
  APPLY_TO_COLS = 2 # parameter for base::apply
  read.data     = readSamples(cov.folder, sample.file, nSamples=5)
  methylDataRaw = read.data$methyldata
  data.env$sampleInfo = read.data$samples
  
  # In each col, find the 90th percentile
  sample.quantiles = apply(totalReads(methylDataRaw), APPLY_TO_COLS, function(col) quantile(col, 1-default.high.exp.filter))
  
  # Function to test if the value in each cell is less than the 90th percentile of its column
  areRowCellsBelowColQuantiles = function(row, col.quantiles){
    all(sapply(1:length(row), function(i) row[i]<col.quantiles[i] )) 
  }
  
  # Create a row predicate on coverage per sample
  valid.counts = apply(totalReads(methylDataRaw), APPLY_TO_ROWS, areRowCellsBelowColQuantiles, sample.quantiles)

  # Create a row predicate on mimumum read count for each sample
  valid.min = apply(totalReads(methylDataRaw), APPLY_TO_ROWS, function(row) all(row>=min.reads ))
  
  # Combine the predicates
  valid.all = valid.counts & valid.min

  # Apply the predicate to the read data
  meth.reads  = methReads(methylDataRaw)[valid.all,]
  total.reads = totalReads(methylDataRaw)[valid.all,]
  cat("Retained", nrow(total.reads),"loci\n")
  
  # Ensure methylated fractions of zero and 1 are never possible
  data.env$b_values = (meth.reads+1)/ (total.reads+2)
  colnames(data.env$b_values) = data.env$sampleInfo$Sample.Name
  data.env$methylDataRaw = methylDataRaw[valid.all,] # keep the valid ranges in the data environent
}

loadData(default.min.reads)


# Mean centre the b-values to prevent fully methylated or unmethylated values
# swamping the results
centred_m = as.data.frame(data.env$b_values) %>% 
  mutate(meanrow = rowMeans(.)) %>% 
  mutate_all(funs(.-meanrow)) %>%
  dplyr::select(-meanrow)
centred_m = as.matrix(centred_m)

# Calculate the cg mean b_values per group
meanLow  = rowMeans(data.env$b_values[,data.env$sampleInfo$BirthweightCat=="low"])
meanHigh = rowMeans(data.env$b_values[,data.env$sampleInfo$BirthweightCat=="high"])
means    = cbind(meanLow, meanHigh, meanLow-meanHigh)

# Create the real and null model matrix for estimation
mod.real = model.matrix(~BirthweightCat+Sex+Season, data=data.env$sampleInfo)
mod.null = model.matrix(~Sex+Season, data=data.env$sampleInfo)

cat("Estimating covariates\n")

# Apply SVA to the data
# Find the number of factors to be estimated, then estimate them
# n.sv  = sva::num.sv(centred_m,mod.real,method="leek")
nsvobj = sva::sva(centred_m,mod.real,mod.null,n.sv=NULL)
if(!exists("nsvobj")){
  cat("SVA failed, estimating number of SVs")
  n.sv  = sva::num.sv(centred_m,mod.real,method="leek")
  cat("Estimated number of surrogate variables is", n.sv)
  quit(save="no")
}

# Add the estimated factors to the model, and fit the new model
cat("\nFitting covariates\n")
mod.real.sv = cbind(mod.real,svobj$sv)
mod.null.sv = cbind(mod.null,svobj$sv)
data.env$fit = limma::lmFit(centred_m, mod.real.sv, method="robust")

# Use ebayes to calculate the test statistics 
cat("Calculating test statistics\n")
fit.e1  = eBayes(data.env$fit)

# Get FDR corrected differential probe list, with any probes below the pvalue threshold
tab.all     = topTable(fit.e1, p.value=0.05, number=Inf, adjust = "fdr")
tab1.birth  = topTable(fit.e1, coef="BirthweightCatlow",p.value=0.05, number=Inf, adjust = "fdr")
tab1.sex    = topTable(fit.e1, coef="SexM",             p.value=0.05, number=Inf, adjust = "fdr")
tab1.season = topTable(fit.e1, coef="Seasonwet",        p.value=0.05, number=Inf, adjust = "fdr")

#' Export a table of significant differentially methylated loci
#' 
#' Adds genomic ranges to the table of loci and exports to the given file.
#' Ignores any loci for which the difference in methylation between groups
#' is less than the given minimum.
#' 
#' @param loci the table of significantly differentially methylated loci
#' @param filename the file to export to
#' @param min_beta_diff the minimum difference between groups
#'
exportSignificantLoci = function(loci, filename, min_beta_diff){
  # Get the names of the significant rows and fetch the corresponding genomic ranges
  sig.rows   = as.numeric(rownames(loci))
  sig.ranges = data.env$methylDataRaw[sig.rows]
  
  # Build the ranges into a table
  result.table = data.frame( chr    = seqnames(sig.ranges),
                             start  = start(sig.ranges),
                             end    = end(sig.ranges),
                             strand = strand(sig.ranges))
  
  read.counts = totalReads(sig.ranges)
  colnames(read.counts) = paste0(data.env$sampleInfo$Sample.Name, "_total_reads")
  
  # Get the proxy beta values for the significant samples
  sig.b_vals = data.env$b_values[sig.rows,]
  colnames(sig.b_vals) = data.env$sampleInfo$Sample.Name
  
  sig.means = means[sig.rows,]
  colnames(sig.means) = c("Mean_Low", "Mean_high", "Diff")
  
  # Mung it all together in one big table
  result.table = cbind(result.table, loci, sig.means, sig.b_vals, read.counts)
  
  # Filter by beta difference
  result.table = result.table %>% filter( abs(Diff)>=min_beta_diff)
  
  min_read_count = min(totalReads(sig.ranges))
  
  table.file = paste0(root.dir, "bms41/Humans/HeroG/sva/", filename, ".", min_read_count, "_reads.", min_beta_diff, "_diff.csv")
  write.table(result.table, file=table.file, sep=",", row.names = F, col.names = T)
}

tables = list("SVA.birth"=tab1.birth, "SVA.sex"=tab1.sex, "SVA.season"=tab1.season)
invisible(mapply(exportSignificantLoci, tables, names(tables), default.min.beta.diff))
