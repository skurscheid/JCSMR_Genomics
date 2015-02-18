#  _TF1_H2A.Bbd_occupancy_data_preparation.R
#
#  Copyright 2015 Sebastian Kurscheid <sebastian.kurscheid@anu.edu.au>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
# load libraries
library("biomaRt")
library("gdata")
library("GO.db")
library("GenomicFeatures")
library("GenomicAlignments")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("rtracklayer")
library("AnnotationDbi")
library("BSgenome.Hsapiens.UCSC.hg19")

# load hg19 based Ensembl data # do this only once, save RData objects than recycle
# human <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
# turns out reads were aligned to hg38...
human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
filters <- listFilters(human)
hsapEnsembl <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
chromInfo <- getChromInfoFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# only use canonical chromosomes, exclude contigs
seqlevels(hsapEnsembl, force = TRUE) <- c(seq(1,22,1), "X", "Y", "MT")

#genome <- BSgenome.Hsapiens.UCSC.hg19
#chromInfo <- seqlengths(genome)[1:25]
#names(chromInfo) <- gsub("chr", "", names(chromInfo))
#names(chromInfo)[25] <- "MT"

save(human, file = "/Volumes/gduserv/Data/Annotations/hsapiens_gene_ensembl_GRCh38.rda")
saveDb(hsapEnsembl, file = "/Volumes/gduserv/Data/Annotations/hsapiens_gene_ensembl_GRCh38_TxDB.sqlite")
# save(chromInfo, file = "/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_chromInfo_GRCh37.rda")

#------------------------------------------------------------------------------
# Prepare GRanges objects for the different genomic regions of interest
#------------------------------------------------------------------------------
gr.tx <- transcripts(hsapEnsembl)

# extract 5' UTRs
grl.5UTR <- fiveUTRsByTranscript(hsapEnsembl)

# reorder GRL per chromosome
sfExport("grl.5UTR")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.5UTR[which(seqnames(grl.5UTR) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.5UTR <- grl

# extract 3' UTRs
grl.3UTR <- threeUTRsByTranscript(hsapEnsembl)

# reorder GRL per chromosome
sfExport("grl.3UTR")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.3UTR[which(seqnames(grl.3UTR) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.3UTR <- grl

# extract CDS
grl.cds <- cdsBy(hsapEnsembl, by = "gene")

# reorder GRL per chromosome
sfExport("grl.cds")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.cds[which(seqnames(grl.cds) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.cds <- grl

# extract genes
gr.genes <- genes(hsapEnsembl)
gr.genes <- sort(gr.genes)

# extract TSSs
gr.tssUp1000Dn50 <- promoters(gr.genes, upstream = 1000, down = 50)

tss_up1000 <- start(gr.genes) - 999
gr.tssUp1000Wholegene <- GRanges(seqnames(gr.genes), IRanges(tss_up1000, end(gr.genes)), strand = strand(gr.genes))
names(gr.tssUp1000Wholegene) <- names(gr.genes)

gr.tssUp1000Dn1000 <- promoters(gr.genes, upstream = 1000, down = 1000)
gr.tssUp125Dn50 <- promoters(gr.genes, upstream = 125, down = 50)

# extract intergenic regions
gr.intergenic <- gaps(gr.genes)
gr.intergenicExcludeTSSUp1000 <- gaps(gr.tss_up1000_wholegene)

# in order to use GRangesList objects to create IRangesList objects for Views of coverage datat,
# the GRangesList objects have to be collapsed back to a chromosome level organisation, i.e. 
# one list GRanges object per chromosome containing all ranges of interest
# it is computationally expensive, therefore we use snowfall
library(snowfall)
sfInit(parallel = TRUE, cpus = 25)
sfLibrary(biovizBase)
sfLibrary(GenomicRanges)
sfExport("canonicalChr")

# Introns
grl.introns <- intronsByTranscript(hsapEnsembl)
sfExport("grl.introns")
n1 <- sfLapply(grl.introns, function(x) {length(x)})
grl.introns <- grl.introns[names(n1[which(n1 > 0)])]

# reorder GRL per chromosome
sfExport("grl.introns")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.introns[which(seqnames(grl.introns) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.introns <- grl

# Exons
grl.exons.gene <- exonsBy(hsapEnsembl, by = "gene") # returns GRangesList
sfExport("grl.exons.gene")

# reorder GRL to per-chromosome
sfExport("grl.exons.gene")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.exons.gene[which(seqnames(grl.exons.gene) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.exons.gene <- grl

# GRangesList with only 1st exon of every gene
grl.exons1st <- GRangesList(sfLapply(grl.exons.gene, function(x) {
  if (length(x) > 0){
    g <- x[1]
    return(g)
  }
}))

# reorder GRL to per-chromosome
sfExport("grl.exons1st")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.exons1st[which(seqnames(grl.exons1st) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.exons1st <- grl

# GRL with remaining exons/gene (if more than 1)
l1 <- sfLapply(grl.exons.gene, function(x) length(x))
sfExport("l1")
grl.exonsRest <- GRangesList(sfLapply(grl.exons.gene[which(l1 >= 2)], function(x) {
  if (length(x) == 2){
    g <- x[2]
    if (class(g) == "GRanges"){
      return(g)
    }
  } else if (length(x) > 2){
    g <- x[2:length(x)]
    if (class(g) == "GRanges"){
      return(g)
    }
  }
}
))

# reorder GRL to per-chromosome
sfExport("grl.exonsRest")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.exonsRest[which(seqnames(grl.exonsRest) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.exonsRest <- grl

# extract starting positions of all introns to generate GRangesList for Exon->Intron boundaries
# first remove objects withouth introns
grl.exonIntron <- GRangesList(sfLapply(grl.introns, function(x) {
  if (length(x) > 0){ 
    d <- data.frame(chr = seqnames(x), start = start(x) - 49, width = 100, strand = as.character(strand(x)))
    gr <- GRanges(as.character(d[,"chr"]), IRanges(start = as.numeric(d[,"start"]), width = as.numeric(d[,"width"])), strand = as.character(d[,"strand"]))
    return(gr)
  }
}
)
)

sfExport("grl.exonIntron")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.exonIntron[which(seqnames(grl.exonIntron) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.exonIntron <- grl

# extract end positions of all introns to generate GRangesList for Intron->Exon boundaries
sfExport("grl.introns")
grl.intronExon <- GRangesList(sfLapply(grl.introns, function(x) {
  if (length(x) > 0){ 
    d <- data.frame(chr = seqnames(x), start = end(x) - 49, width = 100, strand = as.character(strand(x)))
    gr <- GRanges(as.character(d[,"chr"]), IRanges(start = as.numeric(d[,"start"]), width = as.numeric(d[,"width"])), strand = as.character(d[,"strand"]))
    return(gr)
  }
}
)
)

sfExport("grl.intronExon")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.intronExon[which(seqnames(grl.intronExon) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.intronExon <- grl

grl.intronExonFlank25 <- flank(grl.intronExon, width = 25, both = TRUE)
grl.exonIntronFlank25 <- flank(grl.exonIntron, width = 25, both = TRUE)
grl.intronsFlank25 <- flank(grl.introns, width = 25, both = TRUE)

# dest_dir <- "/home/skurscheid/Data/Annotations/hg19"
dest_dir <- "/home/skurscheid/Data/Annotations/hg38"
save(gr.tx, file = paste(dest_dir, "gr.tx.rda" , sep = "/"))
save(grl.exons.gene, file = paste(dest_dir, "grl.exons.gene.rda", sep = "/")) # sorting OK
save(grl.introns, file = paste(dest_dir, "grl.introns.rda", sep = "/")) # sorting OK?
save(grl.5UTR, file = paste(dest_dir, "grl.5UTR.rda", sep = "/")) #
save(grl.3UTR, file = paste(dest_dir, "grl.3UTR.rda", sep = "/")) #
save(grl.cds, file = paste(dest_dir, "grl.cds.rda", sep = "/")) #
save(gr.genes, file = paste(dest_dir, "gr.genes.rda", sep = "/")) # sorting OK?
save(gr.tssUp1000Dn50, file = paste(dest_dir, "gr.tssUp1000Dn50.rda", sep = "/"))
save(gr.tssUp1000Dn1000, file = paste(dest_dir, "gr.tssUp1000Dn1000.rda", sep = "/"))
save(gr.tssUp1000Wholegene, file = paste(dest_dir, "gr.tssUp1000Wholegene.rda", sep = "/"))
save(gr.intergenic, file = paste(dest_dir, "gr.intergenic.rda", sep = "/"))
save(gr.intergenicExcludeTSSUp1000, file = paste(dest_dir, "gr.intergenicExcludeTSSUp1000.rda", sep = "/"))
save(grl.intronExon, file = paste(dest_dir, "grl.intronExon.rda", sep = "/")) # sorting OK
save(grl.intronExonFlank25, file = paste(dest_dir, "grl.introExonFlank25.rda", sep = "/"))
save(grl.exonIntronFlank25, file = paste(dest_dir, "grl.exonIntronFlank25.rda", sep = "/"))
save(grl.intronsFlank25, file = paste(dest_dir, "grl.intronsFlank25.rda", sep = "/"))
save(grl.exons1st, file = paste(dest_dir, "grl.exons1st.rda", sep = "/")) # sorting OK
save(grl.exonsRest, file = paste(dest_dir, "grl.exonsRest.rda", sep = "/")) # sorting OK
save(grl.exonIntron, file = paste(dest_dir, "grl.exonIntron.rda", sep = "/")) # sorting OK


# load the whole bigWig file
data_dir <- "/home/skurscheid/Data/Tremethick/Alignments"
files <- c(paste(data_dir, "/Sample_TF1_1/TF1_1_L001_fill.bw", sep = ""),
           paste(data_dir, "/Sample_TF1_2/TF1_2_combined_L001_fill.bw", sep = ""),
           paste(data_dir, "/Sample_TF1_3/TF1_3_L001_fill.bw", sep = ""),
           paste(data_dir, "/Sample_TF1a_1/TF1a_1_combined_L001_fill.bw", sep = ""),
           paste(data_dir, "/Sample_TF1a_2/TF1a_2_combined_L001_fill.bw", sep = ""),
           paste(data_dir, "/Sample_TF1a_3/TF1a3_combined_L001_fill.bw", sep = ""))

gc()

sapply(files, function(x) print(x))
sapply(files, function(x) file.exists(x))

seq_data <- sapply(files, function(x) import(x))
seq_data_TF1 <- seq_data[1:3]
seq_data_TF1a <- seq_data[4:6]

save(seq_data_TF1, file = "/home/skurscheid/Data/Tremethick/Alignments/seq_data_TF1.rdata")
save(seq_data_TF1a, file = "/home/skurscheid/Data/Tremethick/Alignments/seq_data_TF1a.rdata")

# manually extracted from alignment stats
libSize <- data.frame(sample = c("Sample_TF1_1", "Sample_TF1_2", "Sample_TF1_3", "Sample_TF1a_1", "Sample_TF1a_2", "Sample_TF1a_3"), 
                      size = c(54375559, 49570361, 46122170, 52260493, 45097103, 44502715))
