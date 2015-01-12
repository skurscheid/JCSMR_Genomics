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
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

library("rtracklayer")

# load hg19 based Ensembl data
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters <- listFilters(human)
txdb.hsapEnsembl <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# set up txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# for test purposes only use a small chromosome, e.g. chr12
seqlevels(txdb, force = TRUE) <- c("chr12")

# make GRanges object of transcripts and other features on chr12
gr.tx.chr12 <- transcripts(txdb)
gr.5UTR.chr12 <- fiveUTRsByTranscript(txdb)


# load TF-1 sample 1 from bigWig
# for now only chr12
which <- GRanges(c("12"), IRanges(1,121257530))
gr.tf1.s1 <- import("/Volumes//LaCie//Project_SN877_0258_YWu_JCSMR_human_ChIPseq/Sample_TF1_1/TF1_1_L001_fill.bw", which = which)
