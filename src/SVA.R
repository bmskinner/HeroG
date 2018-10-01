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
  
  source("readMethylationSamples.R") # simplify reading of bismark files
}

setUpWorkSpace()

data.env = new.env() # data for the ongoing analysis

# Variables to change before running on cluster
# root.dir =  "Y:/"
root.dir =  "/mnt/research2/"

sample.file = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra.txt")
cov.folder  = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/cov_files")

data.file = paste0(root.dir, "bms41/Humans/HeroG/data.env.Rdata")

loadData = function(){
  
  read.data     = suppressWarnings(readSamples(cov.folder, sample.file))
  methylDataRaw = read.data$methyldata
  data.env$sampleInfo    = read.data$samples

  # Create the table for SVA input.  The data should be a matrix with features
  # in the rows and samples in the columns
  # Remove rows where the total read count is less than <min.reads> in any sample - not informative

  min.reads = 5
  cat("Filtering reads for", min.reads, "total read depth\n")
  valid_total = apply( totalReads(methylDataRaw), 1, function(row) all(row >=min.reads ))
  meth_reads  = methReads(methylDataRaw)[valid_total,]
  total_reads = totalReads(methylDataRaw)[valid_total,]
  cat("Retained", nrow(total_reads),"reads\n")
  
  data.env$b_values = (meth_reads+1)/ (total_reads+2)
  colnames(data.env$b_values) = data.env$sampleInfo$Sample.Name
  data.env$methylDataRaw = methylDataRaw[valid_total,] # keep the valid ranges in the data environent to be saved
  save(data.env, file = data.file)
  cat("Saved data environment to file\n")
}

if(file.exists(data.file)){
  cat("Loading data environment from file\n")
  load(data.file, envir = globalenv())
} else {
  loadData()
}

# Mean centre the b-values
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
# n.sv  = sva::num.sv(beta_values,mod.real,method="leek")
svobj = sva::sva(centred_m,mod.real,mod.null,n.sv=NULL)

# Add the estimated factors to the model, and fit the new model
cat("\nFitting covariates\n")
mod.real.sv = cbind(mod.real,svobj$sv)
mod.null.sv = cbind(mod.null,svobj$sv)
data.env$fit = limma::lmFit(centred_m, mod.real.sv, method="robust")

save(data.env, file = data.file)
cat("Saved data environment to file\n")

# Use ebayes to calculate the test statistics 
cat("Calculating test statistics\n")
fit.e1  = eBayes(data.env$fit)

# Get FDR corrected differential probe list, with any probes below the pvalue threshold
tab.all     = topTable(fit.e1, p.value=0.05, number=Inf, adjust = "fdr")
tab1.birth  = topTable(fit.e1, coef="BirthweightCatlow",p.value=0.05, number=Inf, adjust = "fdr")
tab1.sex    = topTable(fit.e1, coef="SexM",p.value=0.05, number=Inf, adjust = "fdr")
tab1.season = topTable(fit.e1, coef="Seasonwet",p.value=0.05, number=Inf, adjust = "fdr")

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
  
  table.file = paste0(root.dir, "bms41/Humans/HeroG/", filename, ".", min_read_count, "_reads.", min_beta_diff, "_diff.csv")
  write.table(result.table, file=table.file, sep=",", row.names = F, col.names = T)
}

tables = list(tab1.birth, tab1.sex, tab1.season)
fnames = list("SVA.birth", "SVA.sex", "SVA.season")

invisible(mapply(exportSignificantLoci, tables, fnames, 0.01))




# Estimation by isva

# phenotype = model.matrix(~BirthweightCat, data=data.env$sampleInfo)[,2]
# confounds = model.matrix(~Sex+Season, data=data.env$sampleInfo)[,-1]
# 
# isva.result.nofactors = DoISVA(centred_m, phenotype, pvthCF = 0.01,
#                      th = 0.05, ncomp = NULL,icamethod="JADE")

# isva.result.confounds = DoISVA(centred_m, phenotype, cf.m = confounds, factor.log = c(T, T), pvthCF = 0.01,
#        th = 0.05, ncomp = NULL,icamethod="JADE")

# print(cor(isva.result$isv,confounds))
