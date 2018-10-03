# Read the cgs identified by SVA and compare to BiSeq's DMRs

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
  
  
  biopackages = c("BiSeq", "GenomicRanges", "GenomeGraphs", "biomaRt", "BiSeq", "rtracklayer")
  lapply(biopackages, installIfNeededFromBioconductor)
  
  packages = c("dplyr", "ggplot2", "ggbio")
  lapply(packages, installIfNeeded)
  
  source("readMethylationSamples.R") # simplify reading of bismark files
}
##cq added ggplot
setUpWorkSpace()

# Generate a human reference genome
createReferenceGenome = function(){
  data(UCSC.HG38.Human.CytoBandIdeogram)
  
  # Find the sequence lengths for the chromosomes
  chrlengths = c()
  chrnames   = c()
  for( s in unique(UCSC.HG38.Human.CytoBandIdeogram$Chromosome)){
    chr = UCSC.HG38.Human.CytoBandIdeogram[UCSC.HG38.Human.CytoBandIdeogram$Chromosome == s,]
    end = max(chr$chromEnd)
    cat(s, ": ", end, "\n")
    chrlengths = c(chrlengths, end)
    chrnames = c(chrnames, s)
  }
  hg38.ideo = with(UCSC.HG38.Human.CytoBandIdeogram, GRanges(seqnames = Chromosome, 
                                                             IRanges(start=chromStart+1, 
                                                                     end=chromEnd),
                                                             seqinfo=Seqinfo(chrnames, chrlengths)))
  
  hg38.ideo
}


hg38.ideo = createReferenceGenome()
chrnames = seqlevels(hg38.ideo)
chrlengths = seqlengths(hg38.ideo)
# root.dir =  "/mnt/research2/" # linux
root.dir =  "Y:/" # windows

cov.folder  = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/cov_files")
sample.file = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra.txt")
sampleInfo = readSamples(cov.folder, sample.file, info_only=T)$samples

# Read a list of dmrs from file and build into genomic ranges
dmr.file = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Cord_sample_analysis/DMRs_annotated_kit_male_cord_bloods.txt") 
dmr.data = read.csv(dmr.file, header = T, sep="\t", stringsAsFactors = F)
dmr.data$seqnames = paste0("chr", dmr.data$seqnames)
dmr.granges = with(dmr.data, GRanges(seqnames = seqnames, 
                                     IRanges(start=start, 
                                             end=end),
                                     seqinfo=Seqinfo(seqlevels(hg38.ideo), seqlengths(hg38.ideo))))

#' Import CpGs from file.
#' 
#' Constract a GRanges object with appropriate sequence lengths
#' and drop any sequence names other than main chromsomes
#'
#' @param sva.file the file of CpGs
#'
#' @return the CpGs as a GenomicRanges object
#' @export
#'
createGrangesFromFile = function(sva.file){
  sva.data = read.csv(sva.file, header = T, sep=",", stringsAsFactors = F)
  sva.data$chr = paste0("chr", sva.data$chr)
  
  sva.data = sva.data %>% filter(chr %in% seqnames(hg38.ideo)) # remove unmapped contigs
  
  sva.granges = with(sva.data, GRanges(seqnames = chr, 
                                       IRanges(start=start, 
                                               end=end),
                                       strand = strand,
                                       seqinfo=Seqinfo(seqlevels(hg38.ideo), 
                                                       seqlengths(hg38.ideo))))
  sva.granges
}

#' Find DMRs overlapping CpGs in the given input file
#'
#' @param sva.file the file of CpGs
#' @param dmr.granges a GenomicRanges object containing DMRs
#'
#' @return the CpGs overlapping DMRs as a GenomicRanges object
#' @export
#'
#' @examples
getDmrOverlaps = function(sva.file, dmr.granges){
  
  sva.granges = createGrangesFromFile(sva.file)
  # Allow an offset of 1 in either direction
  sva.granges = GenomicRanges::shift(GenomicRanges::resize(sva.granges, 2), -1)
  

  dmrs.with.sva = subsetByOverlaps(sva.granges,dmr.granges)
  
  
  dmrs.with.sva
}

# Finding CpGs overlapping defined DMRs
birth.dmrs  = getDmrOverlaps(paste0(root.dir, "bms41/Humans/HeroG/SVA.birth.5_reads.0.01_diff.csv"), dmr.granges)
sex.dmrs    = getDmrOverlaps(paste0(root.dir, "bms41/Humans/HeroG/SVA.sex.5_reads.0.01_diff.csv"), dmr.granges)
season.dmrs = getDmrOverlaps(paste0(root.dir, "bms41/Humans/HeroG/SVA.season.5_reads.0.01_diff.csv"), dmr.granges)

# Reading all CpGs
sva.granges    = createGrangesFromFile(paste0(root.dir, "bms41/Humans/HeroG/SVA.birth.5_reads.0.01_diff.csv"))
sva.combined = GenomicRanges::reduce(GenomicRanges::shift(GenomicRanges::resize(sva.granges, 1000),-500))

#' Plot CpGs in samples 
#' 
#' An extra filter is applied: only samples in which the median of one group lies outside the
#' total range of the other group are plotted.
#'
#' @param range.data the GenomicRanges object containing the CpGs
#' @param min_reads the minimim number of reads per sample
#' @param min_diff the minimum methylation difference between high and low birthweight group means
#'
#' @return
#' @export
plotCg = function(range.data, min_reads, min_diff){
  
  chrom = as.character(unique(seqnames(range.data)))
  s = unique(start(range.data))
  e = unique(end(range.data))
  
  filt = sva.data %>% filter(chr==chrom, start>=s, end<=e)

  diff = max(filt$Diff)
  if(diff<min_diff){
    return()
  }

  pval = max(filt$adj.P.Val)
  samples = sampleInfo %>% dplyr::select(Sample.Name, BirthweightCat, Sex, Season)
  values  = as.data.frame(t(filt[,14:49]))
  colnames(values) = c("PctMeth")
  reads   = as.data.frame(t(filt%>%dplyr::select( contains("total_reads"))))
  colnames(reads) = c("TotalReads")
  values$name = gsub("\\.", "-", row.names(values))
  reads$name = gsub("\\.", "-", row.names(reads))
  reads$name = gsub("_total_reads", "", reads$name)
  #reads$fReads = reads$TotalReads / max_reads
  joined = merge(samples, values, by.x = "Sample.Name", by.y="name") %>% group_by(BirthweightCat) %>% mutate(med=median(PctMeth))
  joined = merge(joined, reads, by.x = "Sample.Name", by.y="name")
  
  summ = joined %>% group_by(BirthweightCat) %>% summarise(min=min(PctMeth), max=max(PctMeth), med=median(PctMeth))
  
  m1 = summ[1,]$med
  m2 = summ[2,]$med
  min1 = summ[1,]$min
  min2 = summ[2,]$min
  max1 = summ[1,]$max
  max2 = summ[2,]$max
  
  if( ( m1>min2 && m1<max2 ) && (m2<max1 && m2>min1)  ){
    # skip ranges where the median of the high birthweight is not outside the range of the low birthweight and vice versa
    return()
  }
  cat(chrom, s, "-", e, " : ",m1,"outside (",min2 ,"-",max2, ") ; or", m2,"outside (",min1 ,"-",max1, ")\n")
  
  position = paste0(chrom,"-", s, "-", e)
  condition = paste0(min_reads,"_reads_",min_diff,"_diff")
  ggplot( joined, aes(x=BirthweightCat, y=PctMeth))+
    ggtitle(paste0(min_reads," min total reads; ",signif(diff,digits=3)," diff; p=", signif(pval,digits=3)))+
    geom_hline(aes(yintercept=median(PctMeth)), size=2)+
    geom_violin(draw_quantiles = c(0.5))+ylab("Fraction methylated reads")+
    xlab(position)+
    geom_jitter(aes(col=Sex, size=TotalReads),width=0.1)+ylim(0,1)
  ggsave(file=paste0(root.dir, "bms41/Humans/HeroG/plots/",position,"_", condition,"_values.png"))
}

# match reads and methylation diff ie 10, 0.01 etc
# iterateCg = function(range.data, i){
#   plotCg(range.data[i,], 5, 0.01)
# }


#
#' Find genes overlapping a set of genomic ranges
#' 
#' Ranges with overlapping genes are exported to the given file. 
#'
#' @param ranges.to.annotate the GenomicRanges object to be annotated
#' @param file.name the output file for the annotated ranges. The folder path must exist.
#' @export
#'
#' @examples
findOverlapsWithGenes = function(ranges.to.annotate, file.name){

  output.dir   = paste0(root.dir, "bms41/Humans/HeroG/dmrs/")
  gtf.file     = '/mnt/cgs-fs3/Sequencing/Genome/Human/gtf/ensembl/grch38/release-83/Homo_sapiens.GRCh38.83.gtf'
  gtf.data     = import.gff(gtf.file, format="gtf",feature.type="gene")

  overlaps     = findOverlaps(ranges.to.annotate, gtf.data)
  sva.overlaps = ranges.to.annotate[queryHits(overlaps)]
  gtf.overlaps = gtf.data[subjectHits(overlaps)]
  
  sva.overlaps$annotation = with(gtf.overlaps, paste(gene_id, type, gene_biotype, gene_name, sep='; '))

  cpgs.annotated = as.data.frame(sva.overlaps) %>% 
    group_by(seqnames, start, end, width, strand) %>% 
    summarise(annot = paste(annotation, collapse="|")) 

  write.table(cpgs.annotated, file=file.name, sep="\t", row.names = F)
}

output.dir   = paste0(root.dir, "bms41/Humans/HeroG/dmrs/")
#cpgs overlapping all genes 
#findOverlapsWithGenes(sva.granges, paste0(output.dir, "CpGs_overlapping_annotated_genes_5_reads_example.tsv"))
#overlap with birth dmrs for all cordbloods
# findOverlapsWithGenes(birth.dmrs, paste0(output.dir, "CpGs_overlapping_annotated_genes_5_reads_0.01_example_BWCat_cordbloods_male.tsv"))
#findOverlapsWithGenes(birth.dmrs, paste0(output.dir, "CpGs_overlapping_annotated_genes_5_reads_DMRs_blahblah.tsv"))


seqlevels(dmr.granges) = c(seqlevels(dmr.granges), 21, "Y")
seqlevels(dmr.granges) = paste0("chr", seqlevels(dmr.granges))
seqlengths(dmr.granges) = seqlengths(hg38.ideo)


out.folder = "../plots/circos/"

dmr.outfile = paste0(out.folder, "DMRs from BiSeq.png")

p = ggbio() + circle(hg38.ideo, geom = "ideo", fill = "gray70") +
  circle(dmr.granges, geom = "rect", color = "orange", size=1) +
  circle(sva.combined, geom = "rect", color = "blue", size=1) +
  circle(birth.dmrs, geom = "rect", color = "green", size=1) +
  circle(sex.dmrs, geom = "rect", color = "red", size=1) +
  circle(hg38.ideo, geom = "scale", size = 2)
ggsave(p)
