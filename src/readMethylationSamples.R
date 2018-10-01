library(BiSeq)
library(GenomicRanges)
library(dplyr)
library(tidyr)

#' Title Read Bismark coverage files from the given folder
#'
#' @param cov_folder the folder containing coverage files
#' @param sample_info_file the sample information
#' @param info_only if true, skip the cov files, read only the sample info file and order the samples by cov file name.
#' @param nSamples the number of samples to read. Default (0) reads all samples. Other values (n) return the first n samples.
#'
#' @return a list containing the BSRaw methylation object and the sample info table
#' @export
#'
#' @examples 
#' data = readSamples(folder, file)
#' methylDataRaw = data$methyldata
#' sampleInfo = data$samples
readSamples = function(cov_folder, sample_info_file, info_only=F, nSamples=0 ){
  
  sampleInfo = read.table(sample_info_file, sep="\t", header=T, stringsAsFactors = F)
  
  files  = list.files(path=cov_folder, pattern=".cov", full.names=T)
  files  = grep(paste0(c(as.character(sampleInfo$Sample.Name)), collapse='|'),files, value=TRUE)
  
  if(nSamples>0){
    files = files[1:nSamples]
  }
  
  cat("Files:\n")
  cat(paste0(files, collapse = "\n"))
  cat("\n")
  # Get the file names for processing
  covName_coord = gsub(files, pattern = ".*cov_files/", replacement ="" )
  covName_coord = gsub(covName_coord, pattern = "_R1_001_val_1.fq.gz_bismark_bt2_pe.bismark.cov.gz", replacement ="" )
  #correct the sample order so that it matches the file order of the cov file if needed
  sampleInfo = sampleInfo[match(covName_coord, sampleInfo$Sample.Name),]
  
  if(info_only){
    return(list("samples" = sampleInfo))
  }
  
  methylDataRaw = readBismark(files, sampleInfo)
  return(list("methyldata" = methylDataRaw, "samples" = sampleInfo))
}

# 
# 
#' Title Read the locations of cgs in the hg38 genome build
#'
#' @param hg38.hm450.file the file location
#' @param strip.chr if true, remove "chr" strings from seqnames
#'
#' @return A GRanges object of cg locations in 1-based coordinates.
#' @export
#'
#' @examples
readHg38CgLocations = function(hg38.hm450.file, strip.chr = T){
  hg38.infinium = read.csv(file=hg38.hm450.file, header = T, sep="\t", stringsAsFactors = F)
  row.names(hg38.infinium) = hg38.infinium$probeID
  
  hg38.infinium.subset = hg38.infinium %>% dplyr::filter(!is.na(CpG_chrm))
  
  hg38.ranges = with(hg38.infinium.subset, GRanges(seqnames = CpG_chrm, 
                                                   IRanges(start=CpG_beg, 
                                                           end=CpG_end),
                                                   strand = probe_strand))
  names(hg38.ranges) = hg38.infinium.subset$probeID
  hg38.ranges = GenomicRanges::shift(hg38.ranges, 1)
  if(strip.chr){
    seqlevels(hg38.ranges) = gsub("chr", "", seqlevels(hg38.ranges)) # chr not present in many other samples
  }
  return(hg38.ranges)
}

#' Title Map methyl sequence data to an Illumina cg list
#'
#' @param methyl.data the methylation read data as BSRaw
#' @param sample.info the sample information table
#' @param cg.locations a GRanges object with cg locations to extract
#'
#' @return a list containing the total reads and methylated reads at each cg
#' @export
#'
#' @examples
getReadsOverlappingIlluminaCgs = function(methyl.data, sample.info, cg.locations){
  overlaps = GenomicRanges::findOverlaps(methyl.data, cg.locations, ignore.strand=F)
  
  # Get the overlapping rows
  meth.hits  = methyl.data[queryHits(overlaps)] # select the methyl columns with hits
  illumina.hits = cg.locations[subjectHits(overlaps)]   # select the Illumina cg columns with hits
  
  agggregateReads = function(readMatrix){
    reads = as.data.frame(readMatrix) %>% 
      dplyr::mutate(cgid = names(illumina.hits@ranges)) %>%
      dplyr::group_by(cgid) %>% 
      dplyr::mutate_all(funs( total = sum(.))) %>%
      dplyr::ungroup() %>%
      dplyr::select(cgid, contains("total")) %>%
      dplyr::distinct()
    colnames(reads) = c("cg", sample.info$Sample.Name)
    return(reads)
  }
  total.reads = agggregateReads(totalReads(meth.hits))
  meth.reads = agggregateReads(methReads(meth.hits))
  return(list("total_reads" = total.reads, "meth_reads" = meth.reads))
}