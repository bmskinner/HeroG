# Read a sample Bismark file to compare blood composition CpG locations

library(BiSeq)
library(GenomicRanges)
library(dplyr)
library(tibble)
library(ggplot2)
library(FDb.InfiniumMethylation.hg19) #GRCh37

# source("http://bioconductor.org/biocLite.R")
# biocLite("FlowSorted.CordBloodNorway.450k")
# biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
# library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(bisect)

# root.dir =  "Y:/"
root.dir =  "/mnt/research2/"

sampleInfo = read.table( paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra.txt"), sep="\t", header=T, stringsAsFactors = F)

folder = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/cov_files")
files  = list.files(path=folder, pattern=".cov", full.names=T)
files  = grep(paste0(c(as.character(sampleInfo$Sample.Name)), collapse='|'),files, value=TRUE)

# files = files[1:3]

cat("Files:\n")
cat(paste0(files, collapse = "\n"))
# Get the file names for processing
covName_coord = gsub(files, pattern = ".*cov_files/", replacement ="" )
covName_coord = gsub(covName_coord, pattern = "_R1_001_val_1.fq.gz_bismark_bt2_pe.bismark.cov.gz", replacement ="" )
#correct the sample order so that it matches the file order of the cov file if needed
sampleInfo = sampleInfo[match(covName_coord, sampleInfo$Sample.Name),]

methylDataRaw = readBismark(files, sampleInfo)


cat("Fetching cg ids\n")
# Get the list of cg ids in the bisect reference blood sample
# Ensure it is a list for subsetting later
# Filter the Illumina ranges to just these ids


# Find the hg38 positions of the illumina 450k cgs
readHg38CgLocations = function(){
  hg38.infinium = read.csv( paste0(root.dir, "bms41/Humans/HeroG/hm450.hg38.manifest.tsv/hm450.hg38.manifest"), header = T, sep="\t", stringsAsFactors = F)
  row.names(hg38.infinium) = hg38.infinium$probeID
  
  hg38.infinium.subset = hg38.infinium %>% dplyr::filter(!is.na(CpG_chrm))
  
  hg38.ranges = with(hg38.infinium.subset, GRanges(seqnames = CpG_chrm, 
                                                   IRanges(start=CpG_beg, 
                                                           end=CpG_end),
                                                   strand = probe_strand))
  names(hg38.ranges) = hg38.infinium.subset$probeID
  return(hg38.ranges)
}
hg38.ranges = readHg38CgLocations()

# Get just the cgs that are relevant for adult blood
cg.idlist = unlist(reference_blood[,1])
subset = hg38.ranges[cg.idlist]
seqlevels(subset) = gsub("chr", "", seqlevels(subset)) # chr not present in methyl data

# hg38 uses zero based coordinates
# Bismark cov files uses 1 based coordinates
# Expand by 2bp around the CpG and recentre by shifting
subset = GenomicRanges::shift(subset, 1)


cat("Subsetting\n")
overlaps = findOverlaps(methylDataRaw, subset, ignore.strand=F)

# Get the overlapping rows
meth.hits  = methylDataRaw[queryHits(overlaps)] # select the methyl columns with hits
illumina.hits = subset[subjectHits(overlaps)]   # select the Illumina cg columns with hits
# names(meth.hits@rowRanges) = illumina.hits$names # copy the cg ids into the methyl data

head(totalReads(meth.hits))
head(methReads(meth.hits))


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

write.table(total.reads, file=paste0(root.dir, "bms41/Humans/HeroG/totalReadsAtAdultBloodCgs.csv"), sep=",", row.names = F, col.names = T, quote = F)
write.table(meth.reads, file=paste0(root.dir, "bms41/Humans/HeroG/methReadsAtAdultBloodCgs.csv"), sep=",", row.names = F, col.names = T, quote = F)

checkMaternaContamination = function(){
  
  # The CpGs to detect maternal contamination.
  maternal.cgids = c('cg25556035', 'cg15931839', 'cg02812891', 'cg12634306', 'cg13138089', 
                     'cg15645660', 'cg25241559', 'cg16617301', 'cg19509778', 'cg24767131')
  
  calculatePctMethylated = function(cgid){
    maternal.cgs = infiniumMethylation[cgid]
    maternal.cgs = GenomicRanges::shift(GenomicRanges::resize(maternal.cgs, 4), -1)
    seqlevels(maternal.cgs) = gsub("chr", "", seqlevels(maternal.cgs)) # 'chr' not present in methyl data seqnames
    methOverlaps = IRanges::subsetByOverlaps(methylDataRaw, maternal.cgs, ignore.strand=T)
    cat(cgid, "\n")
    cat(sampleInfo$Sample.Name, collapse="\t", "\n")
    pctMeth = methReads(methOverlaps)/totalReads(methOverlaps)
    cat(totalReads(methOverlaps), "\n")
    cat(methReads(methOverlaps), "\n")
    cat(pctMeth, "\n")
  }
  
  invisible(lapply(maternal.cgids, calculatePctMethylated))
  
  
}

