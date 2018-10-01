# Check for maternal contamination in blood samples
library(BiSeq)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(bisect)
library(tidyr)
library(lazyeval)
source("readMethylationSamples.R") # simplify reading of bismark files

maternal.cgids = c('cg25556035', 'cg15931839', 'cg02812891', 
                   'cg12634306', 'cg13138089', 'cg15645660', 
                   'cg25241559', 'cg16617301', 'cg19509778', 
                   'cg24767131')

maternal.cgid.thresholds = tribble(~cg, ~array.swan, ~array.bmiq, ~pyro,
"cg25556035",  0.123726,  0.085192,  0.057375,
"cg13138089",  0.18849,  0.214598, NA,
"cg12634306",  0.128539,  0.140893, NA,
"cg25241559",  0.209688,  0.230536, NA,
"cg02812891",  0.107235,  0.102181,  0.029425,
"cg15645660",  0.267049,  0.279702, NA,
"cg19509778",  0.117568,  0.062266,NA,
"cg15931839",  0.14165,  0.084144,  0.081525,
"cg16617301",  0.159146,  0.106105, NA,
"cg24767131",  0.145945,  0.185700,NA)

maternal.cgid.thresholds = maternal.cgid.thresholds %>% rowwise() %>% mutate(thresh = min(array.swan, array.bmiq))

root.dir =  "Y:/"
# root.dir =  "/mnt/research2/"

sampleIinfoFile = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra.txt")
cov.folder  = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/cov_files")

read.data     = readSamples(cov.folder, sampleIinfoFile, 3)
methylDataRaw = read.data$methyldata
sampleInfo    = read.data$samples

cat("Fetching cg ids\n")
# Get the list of cg ids in the bisect reference blood sample
# Ensure it is a list for subsetting later
# Filter the Illumina ranges to just these ids


# Find the hg38 positions of the illumina 450k cgs
hg38.ranges = readHg38CgLocations()
subset = hg38.ranges[maternal.cgids]

cat("Subsetting\n")
overlaps = findOverlaps(methylDataRaw, subset, ignore.strand=F)

# Get the overlapping rows
meth.hits  = methylDataRaw[queryHits(overlaps)] # select the methyl columns with hits
illumina.hits = subset[subjectHits(overlaps)]   # select the Illumina cg columns with hits

agggregateReads = function(readMatrix){
  reads = as.data.frame(readMatrix) %>% 
    dplyr::mutate(cgid = names(illumina.hits@ranges)) %>%
    dplyr::group_by(cgid) %>% 
    dplyr::mutate_all(funs( total = sum(.))) %>%
    dplyr::ungroup() %>%
    dplyr::select(cgid, contains("total")) %>%
    dplyr::distinct()
  colnames(reads) = c("cg", sampleInfo$Sample.Name)
  return(reads)
}

total.reads = agggregateReads(totalReads(meth.hits))
meth.reads = agggregateReads(methReads(meth.hits))

write.table(total.reads, file=paste0(root.dir, "bms41/Humans/HeroG/totalReadsAtMaternalCgs.csv"), sep=",", row.names = F, col.names = T, quote = F)
write.table(meth.reads, file=paste0(root.dir, "bms41/Humans/HeroG/methReadsAtMaternalCgs.csv"), sep=",", row.names = F, col.names = T, quote = F)

# read in the table, since the full sample list was run on the cluster

total.reads = read.csv(file=paste0(root.dir, "bms41/Humans/HeroG/totalReadsAtMaternalCgs.csv"), sep=",", header=T, stringsAsFactors = T)
meth.reads = read.csv(file=paste0(root.dir, "bms41/Humans/HeroG/methReadsAtMaternalCgs.csv"), sep=",", header=T, stringsAsFactors = T)
colnames(total.reads) = gsub("\\.", "-", colnames(total.reads))
colnames(meth.reads) = gsub("\\.", "-", colnames(meth.reads))

total.gather = total.reads %>% gather(sampleInfo$Sample.Name, key = "sample_id", value = "total_reads")
meth.gather = meth.reads %>% gather(sampleInfo$Sample.Name, key = "sample_id", value = "meth_reads") %>% select(meth_reads)

combined.gather = cbind(total.gather, meth.gather) %>% mutate(pct.meth = meth_reads/total_reads)

combined.gather.thresh = merge(combined.gather, maternal.cgid.thresholds, by="cg") %>% mutate(is.below.thresh = pct.meth<thresh | is.na(pct.meth))


# Plot the combined data
ggplot(combined.gather.thresh, aes(sample_id, pct.meth))+
  geom_point()+
  geom_hline(data=maternal.cgid.thresholds, aes(yintercept =array.swan))+
  facet_grid(.~cg)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Findsamples with some contamination

test = combined.gather.thresh %>% 
  group_by(sample_id) %>%
  summarise_(n_uncontaminated = interp(~sum(is.below.thresh)))


