# Read the cgs identified by SVA. Find those overlapping genes.
# Export the CpGs and the unique genes to separate files

# Install and load packages
setUpWorkSpace = function(packages, biopackages) {
  
  source("https://bioconductor.org/biocLite.R")
  
  installIfNeededFromBioconductor = function(pkg){
    if(pkg %in% rownames(installed.packages()) == FALSE) {
      biocLite(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
  }
  
  installIfNeeded = function(pkg){
    if(pkg %in% rownames(installed.packages()) == FALSE) {
      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
  }

  lapply(biopackages, installIfNeededFromBioconductor)
  lapply(packages, installIfNeeded)
  
  source("readMethylationSamples.R") # simplify reading of bismark files
}

setUpWorkSpace(packages    = c("dplyr", "ggplot2", "RCircos"),   
               biopackages = c("GenomicRanges", "GenomeGraphs", "ggbio", "rtracklayer" ))

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
working.dir = paste0(root.dir, "bms41/Humans/HeroG/sva/")

gtf.file     = 'O:/Sequencing/Genome/Human/gtf/ensembl/grch38/release-85/Homo_sapiens.GRCh38.85.gtf.gz'
gtf.data     = import.gff(gtf.file, format="gtf",feature.type="gene")
seqlevels(gtf.data) = paste0("chr", seqlevels(gtf.data)) #  add 'chr' prefix to annotation data


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
make.granges.from.sva.file = function(sva.file){
  sva.data = read.csv(sva.file, header = T, sep=",", stringsAsFactors = F)
  sva.data$chr = paste0("chr", sva.data$chr)
  sva.data = sva.data %>% filter(chr %in% seqnames(hg38.ideo)) # remove unmapped contigs
  
  sva.granges = with(sva.data, GRanges(seqnames = chr, 
                                       IRanges(start=start, 
                                               end=end),
                                       strand = strand,
                                       seqinfo=Seqinfo(seqlevels(hg38.ideo), 
                                                       seqlengths(hg38.ideo)),
                                       pval = P.Value))
  sva.granges
}

#' Find genes overlapping a set of genomic ranges
#' 
#' Ranges with overlapping genes are exported to the given file. 
#'
#' @param ranges.to.annotate the GenomicRanges object to be annotated
#' @param file.name the output file for the annotated ranges. The folder path must exist.
#' @export
#'
#' @examples
find.overlapping.genes = function(ranges.to.annotate, file.name){
  
  output.file  = paste0(working.dir, file.name)
  overlaps     = findOverlaps(ranges.to.annotate, gtf.data)
  sva.overlaps = ranges.to.annotate[queryHits(overlaps)]
  gtf.overlaps = gtf.data[subjectHits(overlaps)]
  
  sva.overlaps$annotation = with(gtf.overlaps, paste(gene_id, type, gene_biotype, gene_name, sep='; '))
  
  # Get each CpG annotated with the gene it is within
  cpgs.annotated = as.data.frame(sva.overlaps) %>% 
    group_by(seqnames, start, end, width, strand, pval) %>% 
    summarise(annot = paste(annotation, collapse="|")) 
  
  write.table(cpgs.annotated, file=output.file, sep="\t", row.names = F)
  
  # Get the unique gene list
  genes = as.data.frame(gtf.overlaps) %>% 
    select(seqnames, start, end, gene_id, gene_biotype, gene_name) %>%
    arrange(seqnames, start) %>% distinct()

  write.table(genes, file=paste0(output.file, ".unique.genes.tsv"), sep="\t", row.names = F)
  
}


# Read birth weight associated CpGs
b.file = paste0(working.dir, "SVA.birth.10_reads.0.99_cov.0.01_diff.csv")
b.data = make.granges.from.sva.file(b.file)

# Read sex associated CpGs
s.file = paste0(working.dir, "SVA.sex.10_reads.0.99_cov.0.01_diff.csv")
s.data = make.granges.from.sva.file(s.file)

n.file = paste0(working.dir, "SVA.season.10_reads.0.99_cov.0.01_diff.csv")
n.data = make.granges.from.sva.file(n.file)

# Find overlaps between birthweight and sex
overlaps = findOverlaps(b.data, s.data)
o.data   = b.data[queryHits(overlaps)]

# Subtract sex and season CpGs from birthweight to remove loci 
# resulting from an interaction or either
bw.only.hits = setdiff(b.data, s.data)
bw.only.hits = setdiff(bw.only.hits, n.data)


# Find genes associated with each CpG
find.overlapping.genes(b.data, "CpGs.in.genes.associated.with.birthweight.tsv")
find.overlapping.genes(s.data, "CpGs.in.genes.associated.with.sex.tsv")
find.overlapping.genes(n.data, "CpGs.in.genes.associated.with.season.tsv")
find.overlapping.genes(o.data, "CpGs.in.genes.associated.with.birthweight.and.sex.tsv")
find.overlapping.genes(bw.only.hits, "CpGs.in.genes.associated.with.birthweight.subtracted_sex_season.tsv")


