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
library("GenomeAlignments")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("rtracklayer")

# load hg19 based Ensembl data
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters <- listFilters(human)
hsapEnsembl <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
chromInfo <- getChromInfoFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# for test purposes only use a small chromosome, e.g. chr12
seqlevels(hsapEnsembl, force = TRUE) <- c("12")

# make GRanges object of transcripts and other features on chr12
gr.tx.chr12 <- transcripts(hsapEnsembl)
gr.5UTR.chr12 <- fiveUTRsByTranscript(hsapEnsembl)
gr.exons.chr12 <- exons(hsapEnsembl)

# load TF-1 sample 1 from bigWig
# for now only chr12
which <- GRanges(c("12"), IRanges(1,133851895))
gr.tf1.s1 <- import("/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1_1/TF1_1_L001_fill.bw", which = which)

# use the summarizeOverlaps() function from GenomicAlignments to extract overlapping regions
# sum it up with sum()
sum(data.frame(union =assay(summarizeOverlaps(gr.tx.chr12, gr.tf1.s1)))[,1])
sum(data.frame(union =assay(summarizeOverlaps(gr.exons, gr.tf1.s1)))[,1])
sum(data.frame(union =assay(summarizeOverlaps(gr.5UTR.chr12, gr.tf1.s1)))[,1])
sum(data.frame(union =assay(summarizeOverlaps(which, gr.tf1.s1)))[,1])


# reset TxDB object to include canonical chromosomes
# and create GRanges with regions of interest
seqlevels(hsapEnsembl, force = TRUE) <- c(seq(1,22,1), "X", "Y", "MT")

#------------------------------------------------------------------------------
# Prepare GRanges objects for the different genomic regions of interest
#------------------------------------------------------------------------------
gr.tx <- transcripts(hsapEnsembl)
grl.exons.gene <- exonsBy(hsapEnsembl, by = "gene") # returns GRangesList
grl.introns <- intronsByTranscript(hsapEnsembl)
grl.5UTR <- fiveUTRsByTranscript(hsapEnsembl)
grl.3UTR <- threeUTRsByTranscript(hsapEnsembl)
grl.cds <- cdsBy(hsapEnsembl, by = "gene")
gr.genes <- genes(hsapEnsembl)
gr.genes <- sort(gr.genes)
gr.tss_up1000dn50 <- promoters(gr.genes, upstream = 1000, down = 50)
tss_up1000 <- start(gr.genes) - 1000
gr.tss_up1000_wholegene <- GRanges(seqnames(gr.genes), IRanges(tss_up1000, end(gr.genes)), strand = strand(gr.genes))
names(gr.tss_up1000_wholegene) <- names(gr.genes)
gr.intergenic <- gaps(gr.genes)
gr.intergenic_excludeTSS_up1000 <- gaps(gr.tss_up1000_wholegene)

# extract starting positions of all introns to generate GRangesList for Exon->Intron boundaries
# first remove objects withouth introns
n1 <- lapply(grl.introns, function(x) {length(x)})
grl.introns <- grl.introns[names(n1[which(n1 > 0)])]
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
grl.introns_flank25 <- flank(gr.introns.new, width = 25, both = TRUE)

grl.exons_1st <- GRangesList(lapply(grl.exons.gene, function(x) {
  if (length(x) > 0){
    g <- x[1]
    return(g)
  }
}))

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

# load the whole bigWig file
files <- c("/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1_1/TF1_1_L001_fill.bw",
           "/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1_2/TF1_2_combined_L001_fill.bw",
           "/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1_3/TF1_3_L001_fill.bw",
           "/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1a_1/TF1a_1_combined_L001_fill.bw",
           "/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1a_2/TF1a_2_combined_L001_fill.bw",
           "/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1a_3/TF1a_3_combined_L001_fill.bw")

gr.tf1.s1 <- import("/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1_1/TF1_1_L001_fill.bw")

sapply(files, function(x) print(x))

l1 <- sapply(files, function(x){
  gr.tf1.s1 <- import(x)
  tx <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.tx), gr.tf1.s1)))[,1])
  exons <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exons.gene), gr.tf1.s1)))[,1])
  introns <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.introns), gr.tf1.s1)))[,1])
  UTR5 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.5UTR), gr.tf1.s1)))[,1])
  UTR3 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.3UTR), gr.tf1.s1)))[,1])
  cds <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.cds), gr.tf1.s1)))[,1])
  genes <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.genes), gr.tf1.s1)))[,1])
  tss_up1000dn50 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.tss_up1000dn50), gr.tf1.s1)))[,1]) 
  tss_up1000_wholegene <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.tss_up1000_wholegene), gr.tf1.s1)))[,1]) 
  intergenic <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.intergenic), gr.tf1.s1)))[,1]) 
  intergenic_excludeTSS_up1000 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.intergenic_excludeTSS_up1000), gr.tf1.s1)))[,1]) 
  introns <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.introns), gr.tf1.s1)))[,1]) 
  exon_intron <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exon_intron), gr.tf1.s1)))[,1]) 
  intron_exon <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.intron_exon), gr.tf1.s1)))[,1]) 
  intron_exon_flank25 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.intron_exon_flank25), gr.tf1.s1)))[,1]) 
  exon_intron_flank25 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exon_intron_flank25), gr.tf1.s1)))[,1]) 
  introns_flank25 <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.introns_flank25), gr.tf1.s1)))[,1]) 
  exons_1st <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exons_1st), gr.tf1.s1)))[,1]) 
  exons_rest <- sum(data.frame(union = assay(summarizeOverlaps(reduce(grl.exons_rest), gr.tf1.s1)))[,1]) 
  
  genes_rpb <- genes / sum(unlist(width(gr.genes)))
  tx_rpb <- tx / sum(as.numeric(unlist(width(gr.tx))))
  exons_rpb <- exons / sum(as.numeric(unlist(width(grl.exons.gene))))
  exons_1st_rpb <- exons_1st / sum(as.numeric(unlist(width(grl.exons_1st))))
  exons_rest_rpb <- exons_rest / sum(as.numeric(unlist(width(grl.exons_rest))))
  introns_rpb <- introns / sum(as.numeric(unlist(width(grl.introns))))
  exon_intron_rpb <- exon_intron / sum(as.numeric(unlist(width(grl.exon_intron))))
  intron_exon_flank25_rpb <- intron_exon_flank25 / sum(as.numeric(unlist(width(grl.intron_exon_flank25))))
  exon_intron_flank25_rpb <- exon_intron_flank25 / sum(as.numeric(unlist(width(grl.exon_intron_flank25))))
  introns_flank25_rpb <- introns_flank25 / sum(as.numeric(unlist(width(grl.introns_flank25))))
  UTR5_rpb <- UTR5 / sum(unlist(width(grl.5UTR)))
  UTR3_rpb <- UTR3 / sum(unlist(width(grl.3UTR)))
  cds_rpb <- cds / sum(unlist(width(gr.cds)))
  tss_up1000dn50_rpb <- tss_up1000dn50 / sum(as.numeric(unlist(width(gr.tss_up1000dn50))))
  tss_up1000_wholegene_rpb <- tss_up1000_wholegene / sum(as.numeric(unlist(width(gr.tss_up1000_wholegene))))
  intergenic_rpb <- intergenic / sum(as.numeric(unlist(width(gr.intergenic))))
  intergenic_excludeTSS_up1000_rpb <- intergenic_excludeTSS_up1000 / sum(as.numeric(unlist(width(gr.intergenic_excludeTSS_up1000))))


  df1 <- data.frame(genes_rpb, tx_rpb, exons_rpb, exons_1st_rpb, exons_rest_rpb, introns_rpb, exon_intron_rpb,
                  intron_exon_flank25_rpb, exon_intron_flank25_rpb, introns_flank25_rpb, UTR5_rpb, UTR3_rpb,
                  cds_rpb, tss_up1000dn50_rpb, tss_up1000_wholegene_rpb, intergenic_rpb, intergenic_excludeTSS_up1000_rpb)

  return(df1)
}
)