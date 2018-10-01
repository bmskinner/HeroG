# Read a sample Bismark file to compare blood composition CpG locations

library(BiSeq)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(bisect)
library(tidyr)
root.dir =  "Y:/"
# root.dir =  "/mnt/research2/"

folder = "bms41/Humans/HeroG/"

total.reads = read.table( paste0(root.dir, folder, "totalReadsAtAdultBloodCgs.csv"), sep=",", header=T, stringsAsFactors = F) %>% arrange(cg)
meth.reads  = read.table( paste0(root.dir, folder, "methReadsAtAdultBloodCgs.csv"), sep=",", header=T, stringsAsFactors = F)  %>% arrange(cg)

# sampleInfo = read.table( paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra.txt"), sep="\t", header=T, stringsAsFactors = F)
reference = reference_blood %>% filter(ID %in% total.reads$cg) %>% dplyr::arrange(ID)

# Convert to matrices
m = t(as.matrix(meth.reads[,-1]))
t = t(as.matrix(total.reads[,-1]))
r = as.matrix(reference[,-1])

# Run bisect
results = bisect_supervised(m, t, r, alpha_blood, iterations = 200)

# organizing the results to a data.frame that works with ggplot2
get_visualization_dataframe = function(bisect_results, true_cell_counts) {
  estimates_bin = as.data.frame(bisect_results)
  true_cell_counts = as.data.frame(true_cell_counts)
  
  colnames(estimates_bin) = c("CD4", "CD8", "mono", "Bcells", "NK", "gran")
  colnames(true_cell_counts) = c("CD4", "CD8", "mono", "Bcells", "NK", "gran")
  
  gathered_estimates_bin = estimates_bin %>% gather("CD4", "CD8", "mono", "Bcells", "NK", "gran", key = "cell_type", value = "estimate_norm")
  gathered_truth = true_cell_counts %>% gather("CD4", "CD8", "mono", "Bcells", "NK", "gran", key = "cell_type", value = "truth")
  
  gathered_estimates_bin = gathered_estimates_bin %>% mutate(method = "bin")
  colnames(gathered_estimates_bin) = c("cell_type", "estimate", "method")
  
  estimates = rbind(gathered_estimates_bin)
  truth = rbind(gathered_truth, gathered_truth)
  
  results = cbind(truth, dplyr::select(estimates, "estimate", "method"))
  
  return(results)
}

set.seed(4321)
n_known_samples =36
known_samples_indices = sample.int(nrow(baseline_GSE40279), size = n_known_samples)   
known_samples =  as.matrix(baseline_GSE40279[known_samples_indices, ])

visualization_result = get_visualization_dataframe(results, known_samples)

# plot a scatter plot of true cell types vs estimated.
visualization_result %>% ggplot(aes(truth, estimate, color=cell_type)) + 
  geom_point(size=1, alpha = 0.4) + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("True Cell Proportion") + ylab("Estimated Cell Proportion") + 
  guides(colour = guide_legend(override.aes = list(size=10))) + 
  scale_color_discrete(name = "Cell Type")