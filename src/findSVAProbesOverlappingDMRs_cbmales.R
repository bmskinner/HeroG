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
  
  packages = c("dplyr")
  lapply(packages, installIfNeeded)
  
  source("readMethylationSamples.R") # simplify reading of bismark files
}

setUpWorkSpace()

# root.dir =  "/mnt/research2/"
root.dir =  "Y:/"
# # cq changed below and DMRs to females
cov.folder  = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/cov_files")
sample.file = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/Cord_bloods_low_high_removed_contam_samples_and_outliers_0145n_0397i_with_season_2_plus_extra_variable_males.txt")
sampleInfo = readSamples(cov.folder, sample.file, info_only=T)$samples

dmr.file = paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Cord_sample_analysis/DMRs_annotated_kit_cordblood_males_removed_contam.txt") 
dmr.data = read.csv(dmr.file, header = T, sep="\t", stringsAsFactors = F)
dmr.granges = with(dmr.data, GRanges(seqnames = seqnames, 
                                     IRanges(start=start, 
                                             end=end)))

getDmrOverlaps = function(sva.file){
  sva.data = read.csv(sva.file, header = T, sep=",", stringsAsFactors = F)
  
  sva.granges = with(sva.data, GRanges(seqnames = chr, 
                                       IRanges(start=start, 
                                               end=end),
                                       strand = strand))
  
  # Allow an offset of 1 in either direction
  sva.granges = GenomicRanges::shift(GenomicRanges::resize(sva.granges, 2), -1)
  

  dmrs.with.sva = subsetByOverlaps(sva.granges,dmr.granges)
  dmrs.with.sva
}

birth.dmrs = getDmrOverlaps(paste0(root.dir, "bms41/Humans/HeroG/SVA.birth.5_reads.0.01_diff.csv"))
sex.dmrs = getDmrOverlaps(paste0(root.dir, "bms41/Humans/HeroG/SVA.sex.5_reads.0.01_diff.csv"))
season.dmrs = getDmrOverlaps(paste0(root.dir, "bms41/Humans/HeroG/SVA.season.5_reads.0.01_diff.csv"))


sva.data    = read.csv(paste0(root.dir, "bms41/Humans/HeroG/SVA.birth.5_reads.0.01_diff.csv"), header = T, sep=",", stringsAsFactors = F)
sva.granges = with(sva.data, GRanges(seqnames = chr, 
                                     IRanges(start=start, 
                                             end=end)))

sva.combined = GenomicRanges::reduce(GenomicRanges::shift(GenomicRanges::resize(sva.granges, 1000),-500))

mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org") #connect to BioMart

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
  ##cq samples = sampleInfo %>% dplyr::select(Sample.Name, BirthweightCat, Sex, Season)
  samples = sampleInfo %>% dplyr::select(Sample.Name, BirthweightCat, Season)
  values  = as.data.frame(t(filt[,14:49]))
  colnames(values) = c("PctMeth")
  reads   = as.data.frame(t(filt%>%dplyr::select( contains("total_reads"))))
  colnames(reads) = c("TotalReads")
  values$name = gsub("\\.", "-", row.names(values))
  reads$name = gsub("\\.", "-", row.names(reads))
  reads$name = gsub("_total_reads", "", reads$name)
  reads$fReads = reads$TotalReads / max_reads
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
    folder=paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/plots_cordblood_males_bens_overlap/")
    dir.create(folder)
    ggsave(file=paste0(folder,position,"_", condition,"_values.png"))
  ##ggsave(file=paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/Extra_cord_bloods/plots_cordblood_females_bens_overlap/",position,"_", condition,"_values.png"))
  
  ## cq ggsave(file=paste0(root.dir, "bms41/Humans/HeroG/plots/",position,"_", condition,"_values.png"))
}


#ggsave(file=paste0(root.dir, "crq20/methylSeq/pipelineTest/heroG/yearolds/plots_yearolds_birth_bens_overlap/",position,"_", condition,"_values.png"))
## cq ggsave(file=paste0(root.dir, "bms41/Humans/HeroG/plots/",position,"_", condition,"_values.png"))
iterateCg = function(i){
  plotCg(sva.granges[i,], 5, 0.01)
}

#invisible(lapply(1:length(sva.granges), iterateCg))



findOverlapsWithGenes = function(){

  output.dir   = paste0(root.dir, "bms41/Humans/HeroG/dmrs/")
  ##output.dir   = paste0(root.dir, "bms41/Humans/HeroG/plots/")
  gtf.file     = '/mnt/cgs-fs3/Sequencing/Genome/Human/gtf/ensembl/grch38/release-83/Homo_sapiens.GRCh38.83.gtf'
  gtf.data     = import.gff(gtf.file, format="gtf",feature.type="gene")

  overlaps     = findOverlaps(ranges.to.annotate, gtf.data)
  sva.overlaps = ranges.to.annotate[queryHits(overlaps)]
  gtf.overlaps = gtf.data[subjectHits(overlaps)]
  
  sva.overlaps$annotation = with(gtf.overlaps, paste(gene_id, type, gene_biotype, gene_name, sep='; '))

  cpgs.annotated = as.data.frame(sva.overlaps) %>% 
    group_by(seqnames, start, end, width, strand) %>% 
    summarise(annot = paste(annotation, collapse="|")) 

  #write.table(cpgs.annotated, file=paste0(output.dir, "CpGs_overlapping_annotated_genes_5_reads.tsv"), sep="\t", row.names = F)
  write.table(cpgs.annotated, file=file.name, sep="\t", row.names = F)
}
output.dir   = paste0(root.dir, "bms41/Humans/HeroG/dmrs/")
#ALL GENES
#findOverlapsWithGenes(sva.granges, paste0(output.dir, "CpGs_overlapping_annotated_genes_5_reads_0.01_placenta.tsv"))
#THEN
findOverlapsWithGenes(birth.dmrs, paste0(output.dir, "CpGs_overlapping_annotated_genes_5_reads_0.01_cordbloods_males_BWCat.tsv")) 

  # plotDMR = function(dmr.data){
  #   
  #   temp = subsetByOverlaps(data.env$predictedMeth, dmr.data)
  #   t2   = rowRanges(temp) %>% as.data.frame
  #   t2   = cbind(index=as.numeric(rownames(t2)), t2)
  #   t2$seqnames2 = paste0(t2$seqnames, '_', t2$start)
  #   df = melt(methLevel(temp))
  #   
  #   annot = rowRanges(temp)
  #   annot$seqnames2 = paste0(seqnames(annot), '_', start(annot))
  #   annot = left_join(as.data.frame(annot), DMRs.annotated,by='seqnames2')
  #   annot = paste0(annot[complete.cases(annot),]$annotation, collapse= ' \n ')
  #   
  #   colnames(df) = c('index','Sample.Name','methLevel')
  #   df = left_join(df,as.data.frame(colData(temp)), by='Sample.Name')
  #   df = left_join(df,t2,by='index')
  #   
  #   positions = rowRanges(subsetByOverlaps(data.env$methylDataRaw.filter10.rk.clustered,  dmr.data)) %>% as.data.frame
  #   coverages = totalReads(subsetByOverlaps(data.env$methylDataRaw.filter10.rk.clustered, dmr.data)) %>% as.data.frame
  #   positions = cbind(index=rownames(positions), positions)
  #   coverages = cbind(index=rownames(coverages), coverages)
  #   coverages = melt(coverages)
  #   cov = left_join(coverages,positions,by='index')
  #   cov$seqnames2 = paste0(cov$seqnames, '_', cov$start)
  #   colnames(cov)[2] = 'Sample.Name'
  #   cov = cov %>% select(Sample.Name,value,seqnames2)
  #   
  #   df = left_join(df,cov,by=c('Sample.Name','seqnames2'))
  #   df$seqnames2 = factor(df$seqnames2,levels=unique(df$seqnames2))
  #   
  #   p = ggplot(df,aes(x=seqnames2, y=methLevel, col=value)) +
  #     geom_point() + 
  #     facet_wrap(~Group) + 
  #     ggtitle(annot) + 
  #     scale_fill_gradient() + 
  #     theme(axis.text.x=element_text(angle=90,vjust=0.5))
  #   
  #   #change those 2 lines to change the location where the graphs are saved.
  #   ggsave(plot=p, filename=paste0(outputDir, 'DMR_',df$seqnames2[1],'.png'))
  # }
  # 
  # invisible(lapply(dmrs, plotDMR))
}

# Plotting genomic ranges

# plotCg = function(row.data){
#   
#   cat(row.data$chr, row.data$start, row.data$Diff, "\n")
#   plusStrand = makeGeneRegion(chromosome = row.data$chr, start = row.data$start-1000, end = row.data$end+1000, strand = "+", biomart = mart, dp = DisplayPars(plotId = F, idRotation = 90, idColor = "green"))
#   minStrand = makeGeneRegion( chromosome = row.data$chr, start = row.data$start-1000, end = row.data$end+1000, strand = "-", biomart = mart, dp = DisplayPars(plotId = F, idRotation = 90, idColor = "green"))
#   genomeAxis = makeGenomeAxis(add53=TRUE, add35=TRUE)
# 
#   samples = sampleInfo %>% dplyr::select(Sample.Name, BirthweightCat)
#   values  = as.data.frame(t(row.data[,14:49]))
#   values$name = gsub("\\.", "-", row.names(values))
#   
#   
#   joined = merge(samples, values, by.x = "Sample.Name", by.y="name")
#   
#   allMatrix = data.matrix( t(joined[,3]), rownames.force = NA)
#   matColours = data.matrix( t( ifelse(joined[,2]=="low", "blue", "red")), rownames.force = NA)
# 
#   allArray = makeGenericArray(intensity = allMatrix , probeStart = row.data$start, dp = DisplayPars(color=matColours, type="point"))
# 
#   png(file=paste0(root.dir, "bms41/Humans/HeroG/plots/",row.data$chr,"-",row.data$start,".png"))
#   gdPlot(list(allArray,plusStrand, genomeAxis, minStrand ))
#   dev.off()
# }

# iterateCg = function(i){
#   d = sva.data[i,]
#   if(abs(d$Diff)>0.20){
#     plotCg(d)
#   }
# }


# invisible(lapply(1:nrow(sva.data), iterateCg))

