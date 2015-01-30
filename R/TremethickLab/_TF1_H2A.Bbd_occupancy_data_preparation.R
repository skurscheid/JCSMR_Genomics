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
human <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
filters <- listFilters(human)
hsapEnsembl <- makeTranscriptDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", dataset = "hsapiens_gene_ensembl")
seqlevels(hsapEnsembl, force = TRUE) <- c(seq(1,22,1), "X", "Y", "MT")

genome <- BSgenome.Hsapiens.UCSC.hg19
chromInfo <- seqlengths(genome)[1:25]
names(chromInfo) <- gsub("chr", "", names(chromInfo))
names(chromInfo)[25] <- "MT"

save(human, file = "/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_GRCh37.rda")
saveDb(hsapEnsembl, file = "/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_GRCh37_TxDB.sqlite")
save(chromInfo, file = "/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_chromInfo_GRCh37.rda")

#------------------------------------------------------------------------------
# Prepare GRanges objects for the different genomic regions of interest
#------------------------------------------------------------------------------
gr.tx <- transcripts(hsapEnsembl)
seqlengths(gr.tx) <- chromInfo
grl.exons.gene <- exonsBy(hsapEnsembl, by = "gene") # returns GRangesList
seqlengths(grl.exons.gene) <- chromInfo
grl.introns <- intronsByTranscript(hsapEnsembl)
seqlengths(grl.introns) <- chromInfo
n1 <- lapply(grl.introns, function(x) {length(x)})
grl.introns <- grl.introns[names(n1[which(n1 > 0)])]
grl.5UTR <- fiveUTRsByTranscript(hsapEnsembl)
seqlengths(grl.5UTR) <- chromInfo
grl.3UTR <- threeUTRsByTranscript(hsapEnsembl)
seqlengths(grl.3UTR) <- chromInfo
grl.cds <- cdsBy(hsapEnsembl, by = "gene")
seqlengths(grl.cds) <- chromInfo
gr.genes <- genes(hsapEnsembl)
seqlengths(gr.genes) <- chromInfo
gr.genes <- sort(gr.genes)
gr.tss_up1000dn50 <- promoters(gr.genes, upstream = 1000, down = 50)
tss_up1000 <- start(gr.genes) - 1000
gr.tss_up1000_wholegene <- GRanges(seqnames(gr.genes), IRanges(tss_up1000, end(gr.genes)), strand = strand(gr.genes))
gr.tss_up1000_dn1000 <- GRanges(seqnames(gr.genes), IRanges(tss_up1000, tss_up1000 + 2000), strand = strand(gr.genes))
names(gr.tss_up1000_wholegene) <- names(gr.genes)
names(gr.tss_up1000_dn1000) <- names(gr.genes)
gr.intergenic <- gaps(gr.genes)
gr.intergenic_excludeTSS_up1000 <- gaps(gr.tss_up1000_wholegene)

# extract starting positions of all introns to generate GRangesList for Exon->Intron boundaries
# first remove objects withouth introns
grl.exon_intron <- GRangesList(lapply(grl.introns, function(x) {
  if (length(x) > 0){ 
    d <- data.frame(chr = seqnames(x), start = start(x), width = 1, strand = as.character(strand(x)))
    gr <- GRanges(as.character(d[,"chr"]), IRanges(start = as.numeric(d[,"start"]), width = as.numeric(d[,"width"])), strand = as.character(d[,"strand"]))
    return(gr)
  }
}
)
)

# extract end positions of all introns to generate GRangesList for Intron->Exon boundaries
grl.intron_exon <- GRangesList(lapply(grl.introns, function(x) {
  if (length(x) > 0){ 
    d <- data.frame(chr = seqnames(x), start = end(x), width = 1, strand = as.character(strand(x)))
    gr <- GRanges(as.character(d[,"chr"]), IRanges(start = as.numeric(d[,"start"]), width = as.numeric(d[,"width"])), strand = as.character(d[,"strand"]))
    return(gr)
  }
}
)
)

grl.intron_exon_flank25 <- flank(grl.intron_exon, width = 25, both = TRUE)
grl.exon_intron_flank25 <- flank(grl.exon_intron, width = 25, both = TRUE)
grl.introns_flank25 <- flank(grl.introns, width = 25, both = TRUE)

grl.exons_1st <- GRangesList(lapply(grl.exons.gene, function(x) {
  if (length(x) > 0){
    g <- x[1]
    return(g)
  }
}))

l1 <- lapply(grl.exons.gene, function(x) length(x))
grl.exons_rest <- GRangesList(lapply(grl.exons.gene[which(l1 >= 2)], function(x) {
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

dest_dir <- "/home/skurscheid/Data/Annotations/hg19"
save(gr.tx, file = paste(dest_dir, "gr.tx" , sep = "/"))
save(grl.exons.gene, file = paste(dest_dir, "grl.exons.gene", sep = "/"))
save(grl.introns, file = paste(dest_dir, "grl.introns", sep = "/"))
save(grl.5UTR, file = paste(dest_dir, "grl.5UTR", sep = "/"))
save(grl.3UTR, file = paste(dest_dir, "grl.3UTR", sep = "/"))
save(grl.cds, file = paste(dest_dir, "grl.cds", sep = "/"))
save(gr.genes, file = paste(dest_dir, "gr.genes", sep = "/"))
save(gr.tss_up1000dn50, file = paste(dest_dir, "gr.tss_up1000dn50", sep = "/"))
save(gr.tss_up1000_dn1000, file = paste(dest_dir, "gr.tss_up1000_dn1000", sep = "/"))
save(gr.tss_up1000_wholegene, file = paste(dest_dir, "gr.tss_up1000_wholegene", sep = "/"))
save(gr.intergenic, file = paste(dest_dir, "gr.intergenic", sep = "/"))
save(gr.intergenic_excludeTSS_up1000, file = paste(dest_dir, "gr.intergenic_excludeTSS_up1000", sep = "/"))
save(grl.intron_exon, file = paste(dest_dir, "grl.intron_exon", sep = "/"))
save(grl.intron_exon_flank25, file = paste(dest_dir, "grl.intron_exon_flank25", sep = "/"))
save(grl.exon_intron_flank25, file = paste(dest_dir, "grl.exon_intron_flank25", sep = "/"))
save(grl.introns_flank25, file = paste(dest_dir, "grl.introns_flank25", sep = "/"))
save(grl.exons_1st, file = paste(dest_dir, "grl.exons_1st", sep = "/"))
save(grl.exons_rest, file = paste(dest_dir, "grl.exons_rest", sep = "/"))
save(grl.exon_intron, file = paste(dest_dir, "grl.exon_intron", sep = "/"))


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
