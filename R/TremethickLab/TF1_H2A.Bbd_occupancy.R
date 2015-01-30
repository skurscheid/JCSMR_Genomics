#  <>
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

# load hg19 based Ensembl data # do this only once, save RData objects than recycle
load("/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_GRCh37.rda")
hsapEnsembl <- loadDb("/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_GRCh37_TxDB.sqlite")
load("/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_chromInfo_GRCh37.rda")

# reset TxDB object to include canonical chromosomes
# and create GRanges with regions of interest
seqlevels(hsapEnsembl, force = TRUE) <- c(seq(1,22,1), "X", "Y", "MT")

# load GRanges objects
dest_dir <- "/home/skurscheid/Data/Annotations/hg19"
load(paste(dest_dir, "gr.tx" , sep = "/"))
load(paste(dest_dir, "grl.exons.gene", sep = "/"))
load(paste(dest_dir, "grl.introns", sep = "/"))
load(paste(dest_dir, "grl.5UTR", sep = "/"))
load(paste(dest_dir, "grl.3UTR", sep = "/"))
load(paste(dest_dir, "grl.cds", sep = "/"))
load(paste(dest_dir, "gr.genes", sep = "/"))
load(paste(dest_dir, "gr.tss_up1000dn50", sep = "/"))
load(paste(dest_dir, "gr.tss_up1000_dn1000", sep = "/"))
load(paste(dest_dir, "gr.tss_up1000_wholegene", sep = "/"))
load(paste(dest_dir, "gr.intergenic", sep = "/"))
load(paste(dest_dir, "gr.intergenic_excludeTSS_up1000", sep = "/"))
load(paste(dest_dir, "grl.intron_exon", sep = "/"))
load(paste(dest_dir, "grl.intron_exon_flank25", sep = "/"))
load(paste(dest_dir, "grl.exon_intron_flank25", sep = "/"))
load(paste(dest_dir, "grl.introns_flank25", sep = "/"))
load(paste(dest_dir, "grl.exons_1st", sep = "/"))
load(paste(dest_dir, "grl.exons_rest", sep = "/"))
load(paste(dest_dir, "grl.exon_intron", sep = "/"))

# load sequencing data (from bigWig files)
load("/home/skurscheid/Data/Tremethick/Alignments/seq_data.rdata")

# manually extracted from alignment stats
libSize <- data.frame(sample = c("Sample_TF1_1", "Sample_TF1_2", "Sample_TF1_3", "Sample_TF1a_1", "Sample_TF1a_2", "Sample_TF1a_3"), 
                      size = c(54375559, 49570361, 46122170, 52260493, 45097103, 44502715))

.summarize <- function(y){
  x <- seq_data[[y]]
  lib <- libSize[y,]$size
  
  tx <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.tx), x)))[,1]) / lib * 1000000
  exons <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exons.gene), x)))[,1]) / lib * 1000000
  introns <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.introns), x)))[,1]) / lib * 1000000
  UTR5 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.5UTR), x)))[,1]) / lib * 1000000
  UTR3 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.3UTR), x)))[,1]) / lib * 1000000
  cds <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.cds), x)))[,1]) / lib * 1000000
  genes <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.genes), x)))[,1]) / lib * 1000000
  tss_up1000dn50 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.tss_up1000dn50), x)))[,1]) / lib * 1000000
  tss_up1000_wholegene <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.tss_up1000_wholegene), x)))[,1]) / lib * 1000000
  intergenic <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.intergenic), x)))[,1]) / lib * 1000000
  intergenic_excludeTSS_up1000 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.intergenic_excludeTSS_up1000), x)))[,1]) / lib * 1000000
  introns <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.introns), x)))[,1]) / lib * 1000000
  exon_intron <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exon_intron), x)))[,1]) / lib * 1000000
  intron_exon <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.intron_exon), x)))[,1]) / lib * 1000000
  intron_exon_flank25 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.intron_exon_flank25), x)))[,1]) / lib * 1000000 
  exon_intron_flank25 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exon_intron_flank25), x)))[,1]) / lib * 1000000
  introns_flank25 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.introns_flank25), x)))[,1]) / lib * 1000000
  exons_1st <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exons_1st), x)))[,1]) / lib * 1000000
  exons_rest <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exons_rest), x)))[,1]) / lib * 1000000
  
  tx_rpb <- tx / sum(as.numeric(unlist(width(gr.tx))))
  exons_rpb <- exons / sum(as.numeric(unlist(width(grl.exons.gene))))
  introns_rpb <- introns / sum(as.numeric(unlist(width(grl.introns))))
  UTR5_rpb <- UTR5 / sum(as.numeric(unlist(width(grl.5UTR))))
  UTR3_rpb <- UTR3 / sum(as.numeric(unlist(width(grl.3UTR))))
  cds_rpb <- cds / sum(as.numeric(unlist(width(grl.cds))))
  genes_rpb <- genes / sum(as.numeric(unlist(width(gr.genes))))
  tss_up1000dn50_rpb <- tss_up1000dn50 / sum(as.numeric(unlist(width(gr.tss_up1000dn50))))
  tss_up1000_wholegene_rpb <- tss_up1000_wholegene / sum(as.numeric(unlist(width(gr.tss_up1000_wholegene))))
  intergenic_rpb <- intergenic / sum(as.numeric(unlist(width(gr.intergenic))))
  intergenic_excludeTSS_up1000_rpb <- intergenic_excludeTSS_up1000 / sum(as.numeric(unlist(width(gr.intergenic_excludeTSS_up1000))))
  introns_rpb <- introns / sum(as.numeric(unlist(width(grl.introns))))
  exon_intron_rpb <- exon_intron / sum(as.numeric(unlist(width(grl.exon_intron))))
  intron_exon_rpb <- intron_exon / sum(as.numeric(unlist(width(grl.intron_exon))))
  intron_exon_flank25_rpb <- intron_exon_flank25 / sum(as.numeric(unlist(width(grl.intron_exon_flank25))))
  exon_intron_flank25_rpb <- exon_intron_flank25 / sum(as.numeric(unlist(width(grl.exon_intron_flank25))))
  introns_flank25_rpb <- introns_flank25 / sum(as.numeric(unlist(width(grl.introns_flank25))))
  exons_1st_rpb <- exons_1st / sum(as.numeric(unlist(width(grl.exons_1st))))
  exons_rest_rpb <- exons_rest / sum(as.numeric(unlist(width(grl.exons_rest))))
  
  df1 <- data.frame(tx_rpb,
                    exons_rpb,
                    introns_rpb,
                    UTR5_rpb,
                    UTR3_rpb,
                    cds_rpb,
                    genes_rpb,
                    tss_up1000dn50_rpb,
                    tss_up1000_wholegene_rpb,
                    intergenic_rpb,
                    intergenic_excludeTSS_up1000_rpb,
                    introns_rpb,
                    exon_intron_rpb,
                    intron_exon_rpb,
                    intron_exon_flank25_rpb,
                    exon_intron_flank25_rpb,
                    introns_flank25_rpb,
                    exons_1st_rpb,
                    exons_rest_rpb)
                      
  return(df1)  
}

require(snowfall)
sfInit(parallel = TRUE, cpus = 6)
sfLibrary(GenomicRanges)
sfLibrary(GenomicAlignments)
sfExport(".summarize")
sfExport("libSize")
sfExport("seq_data")
sfExport("gr.tx")
sfExport("grl.exons.gene")
sfExport("grl.introns")
sfExport("grl.5UTR")
sfExport("grl.3UTR")
sfExport("grl.cds")
sfExport("gr.genes")
sfExport("gr.tss_up1000dn50")
sfExport("gr.tss_up1000_wholegene")
sfExport("gr.intergenic")
sfExport("gr.intergenic_excludeTSS_up1000")
sfExport("grl.introns")
sfExport("grl.exon_intron")
sfExport("grl.intron_exon")
sfExport("grl.intron_exon_flank25")
sfExport("grl.exon_intron_flank25")
sfExport("grl.introns_flank25")
sfExport("grl.exons_1st")
sfExport("grl.exons_rest")

l1 <- sfLapply(seq_along(seq_data), function(x) .summarize(x))
sfStop()

# preparing plotting of barplots with ggplot2
library(reshape)
library(ggplot2)

# create data frame from list
df1 <- rbind(l1[[1]], l1[[2]], l1[[3]], l1[[4]], l1[[5]], l1[[6]])
rownames(df1) <- c("Sample_TF1_1", "Sample_TF1_2", "Sample_TF1_3", "Sample_TF1a_1", "Sample_TF1a_2", "Sample_TF1a_3")

# normalize to intergenic rpm/b
df1.norm <- df1 / df1$intergenic_rpb
df1.norm$id <- rownames(df1.norm)
df1.norm$group <- unlist(lapply(strsplit(df1.norm$id, "\\_"), function(x) x[2]))
df1.norm.melt <- melt.data.frame(df1.norm, id.vars = c("group", "id"))

dfc <- summarySE(df1.norm.melt, measurevar="value", groupvars=c("group", "variable"))

dfc2 <- dfc
dfc2$variable <- factor(dfc2$variable)

p <- ggplot(dfc2, aes(x=variable, y=value, fill=group)) + 
     geom_bar(position=position_dodge(), stat="identity") +
     geom_errorbar(aes(ymin=value-se, ymax=value+se),
                   width=.2,                    # Width of the error bars
                   position=position_dodge(.9))

p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12))
