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
  suppressPackageStartupMessages(source("readMethylationSamples.R")) # simplify reading of bismark files
}

setUpWorkSpace()

data.env = new.env() # data for the ongoing analysis

# Variables to change before running on cluster
# root.dir =  "Y:/"
root.dir =  "/mnt/research2/"

sample.file = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra.txt")
cov.folder  = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/cov_files")

default.min.beta.diff   = 0.01  # minimum required difference in methylation fraction
default.min.reads       = 5     # minimum number of reads in every sample
default.high.cov.filter = 0.99  # exclude loci with more than this proportion of the maximum read depth in any sample
default.model.formula   = formula(~BirthweightCat+Sex+Season)
default.null.formaula   = formula(~Sex+Season)


#' Load and filter data
#' 
#' Create the table for SVA input.  The data should be a matrix with features
#' in the rows and samples in the columns.
#' 
#' Filter the incoming data by readcount. Discard the CpGs with the top n% of total read coverage in each sample; this should
#' help remove artefacts due to PCR duplication
#' 
#' @param  min.reads Only include CpGs for which there are this many reads in every sample
#' @param  cov.filter Exclude loci with more than the given percentile of reads in any sample
#' @return a matrix with CpGs in rows and samples in columns.
loadData = function(min.reads, cov.filter){
  cat("Excluding loci with more than", cov.filter,"read coverage in any sample\n")
  cat("Excluding loci with less than", min.reads, "reads in any sample\n")
  
  APPLY_TO_ROWS = 1 # parameter for base::apply
  APPLY_TO_COLS = 2 # parameter for base::apply
  read.data     = readSamples(cov.folder, sample.file) #, nSamples=3
  cat("Read samples\n")
  methylDataRaw = read.data$methyldata
  data.env$sampleInfo = read.data$samples
  
  starting.reads = totalReads(methylDataRaw)
  names(starting.reads) = data.env$sampleInfo$Sample.Name
  
  # Plot the read counts for the given loci
  plotReadCounts = function(reads, file.name){
    tryCatch({
      cat("Creating coverage plot\n")
      
      df = as.data.frame(reads)
      names(df) = data.env$sampleInfo$Sample.Name
      
      gathered = df %>% gather(data.env$sampleInfo$Sample.Name, key = "sample_id", value = "total_reads")
      plot.file  = paste0(root.dir, "bms41/Humans/HeroG/plots/sva/", file.name ,".png")
      g = ggplot(gathered, aes(x=sample_id, y=total_reads))+
        geom_violin()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      cat("Saving coverage image to", plot.file, "\n")
      ggsave(file=plot.file, plot=g)
    }, error=function(e){
      cat("Error making coverage image\n", paste0(e, collapse = "\n"))
    }, warning=function(e){
      cat("Issue making coverage image\n", paste0(e, collapse = "\n"))
    })
  }
  
  plotReadCounts(starting.reads, "Total reads")
  
  # In each col, find the nth percentile.
  sample.quantiles = apply(starting.reads, APPLY_TO_COLS, function(col) quantile(col, cov.filter, type=8))
  cat("Sample coverage filters:\n", paste0(sample.quantiles, collapse = "\n"), "\n")

  cat("Creating predicate...\n")

  # Predicate for min reads per sample
  valid.min =  starting.reads[,1]>=min.reads
  for(i in 1:ncol(starting.reads)) { valid.min = valid.min & starting.reads[,i]>=min.reads }

  # Predicate for max reads per sample
  valid.max = starting.reads[,1]<sample.quantiles[1]
  for(i in 1:ncol(starting.reads)) { valid.max = valid.max & starting.reads[,i]<sample.quantiles[i] }

  # Combine predicates
  valid.all = valid.min & valid.max

  # Apply the predicate to the read data
  cat("Applying predicate...\n")
  meth.reads  = methReads(methylDataRaw)[valid.all,TRUE]
  total.reads = totalReads(methylDataRaw)[valid.all,TRUE]
  cat("Retained", nrow(total.reads),"of", nrow(starting.reads),"loci\n")
  
  plotReadCounts(total.reads, paste0("Filtered reads min-reads_",min.reads, " max-cov_", cov.filter))
  
  # Ensure methylated fractions of zero and 1 are never possible
  data.env$b_values = (meth.reads+1)/ (total.reads+2)
  colnames(data.env$b_values) = data.env$sampleInfo$Sample.Name
  data.env$methylDataRaw = methylDataRaw[valid.all,TRUE] # keep the valid ranges in the data environent
}

loadData(default.min.reads, default.high.cov.filter)


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
mod.real = model.matrix(default.model.formula, data=data.env$sampleInfo) 
mod.null = model.matrix(default.null.formaula, data=data.env$sampleInfo)

cat("Estimating surrogate variables\n")

# Apply SVA to the data
svobj = sva::sva(centred_m,mod.real,mod.null,n.sv=NULL)
if(!exists("svobj")){
  cat("SVA failed using BE method, estimating number of SVs using Leek method\n")
  n.sv  = sva::num.sv(centred_m,mod.real,method="leek")
  cat("Estimated number of surrogate variables is", n.sv,"\n")
  if(n.sv>0){
    svobj = sva::sva(centred_m,mod.real,mod.null,n.sv=n.sv)
  } else {
    cat("No surrogate variables detected\n")
  }
}

# Add the estimated factors to the model, and fit the new model
cat("\nFitting ",svobj$sv,"surrogate variables\n")
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
    
  table.file = paste0(root.dir, "bms41/Humans/HeroG/sva/", filename, ".", default.min.reads, "_reads.", default.high.cov.filter, "_cov.", min_beta_diff, "_diff.csv")
  write.table(result.table, file=table.file, sep=",", row.names = F, col.names = T)
}

tables = list("SVA.birth"=tab1.birth, "SVA.sex"=tab1.sex, "SVA.season"=tab1.season)


invisible(mapply(exportSignificantLoci, tables, names(tables), default.min.beta.diff))

