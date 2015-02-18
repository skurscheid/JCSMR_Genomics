#  TF1_H2A.Bbd_metagene_plots.R
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
library("GO.db")
library("GenomicFeatures")
library("GenomicAlignments")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("rtracklayer")
library("AnnotationDbi")
library("BSgenome.Hsapiens.UCSC.hg19")
library("metagene")
library("GenomicRanges")
library("RSQLite")

# load hg38 based Ensembl data # do this only once, save RData objects than recycle
load("/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_GRCh38.rda")
hsapEnsembl <- loadDb("/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_GRCh38_TxDB.sqlite")

# reset TxDB object to include canonical chromosomes
# and create GRanges with regions of interest
seqlevels(hsapEnsembl, force = TRUE) <- c(seq(1,22,1), "X", "Y", "MT")

# load GRanges objects containing regions-of-interest
dest_dir <- "/home/skurscheid/Data/Annotations/hg38"
load(paste(dest_dir, "gr.tx.rda" , sep = "/"))
load(paste(dest_dir, "grl.exons.gene.rda", sep = "/"))
load(paste(dest_dir, "grl.introns.rda", sep = "/"))
load(paste(dest_dir, "grl.5UTR.rda", sep = "/"))
load(paste(dest_dir, "grl.3UTR.rda", sep = "/"))
load(paste(dest_dir, "grl.cds.rda", sep = "/"))
load(paste(dest_dir, "gr.genes.rda", sep = "/"))
load(paste(dest_dir, "gr.tssUp1000Dn50.rda", sep = "/"))
load(paste(dest_dir, "gr.tssUp1000Dn1000.rda", sep = "/"))
load(paste(dest_dir, "gr.tssUp1000Wholegene.rda", sep = "/"))
load(paste(dest_dir, "gr.intergenic.rda", sep = "/"))
load(paste(dest_dir, "gr.intergenicExcludeTSSUp1000.rda", sep = "/"))
load(paste(dest_dir, "grl.intronExon.rda", sep = "/"))
load(paste(dest_dir, "grl.introExonFlank25.rda", sep = "/"))
load(paste(dest_dir, "grl.exonIntronFlank25.rda", sep = "/"))
load(paste(dest_dir, "grl.intronsFlank25.rda", sep = "/"))
load(paste(dest_dir, "grl.exons1st.rda", sep = "/"))
load(paste(dest_dir, "grl.exonsRest.rda", sep = "/"))
load(paste(dest_dir, "grl.exonIntron.rda", sep = "/"))

# retain only canonical chromosomes, exclude other assembly artefacts
canonicalChr <- c(seq(1,22,1), "X", "Y", "MT")

# prepare GRangesList objects
grl.tssUp125Dn50 <- GRangesList(sapply(canonicalChr, function(x) gr.tssUp125Dn50[which(seqnames(gr.tssUp125Dn50) %in% x)]))
grl.tssUp1000Dn1000 <- GRangesList(sapply(canonicalChr, function(x) gr.tssUp1000Dn1000[which(seqnames(gr.tssUp1000Dn1000) %in% x)]))

#-----------------------------------------------------------------------------------
# BAM file processing - need to transfer into separate script
#-----------------------------------------------------------------------------------
bamDir <- "/home/skurscheid/Data/Tremethick/Alignments/"
bamFiles <- c(paste(bamDir, "Sample_TF1_1/TF1_1_L001_fill_coordinate_sorted_q20.bam", sep = ""),
              paste(bamDir, "Sample_TF1_2/TF1_2_combined_L001_fill_coordinate_sorted_q20.bam", sep = ""),
              paste(bamDir, "Sample_TF1_3/TF1_3_L001_fill_position_sorted_q20.bam", sep = ""),
              paste(bamDir, "Sample_TF1a_1/TF1a_1_combined_L001_fill_coordinate_sorted_q20.bam", sep = ""),
              paste(bamDir, "Sample_TF1a_2/TF1a_2_combined_L001_fill_coordinate_sorted_q20.bam", sep = ""),
              paste(bamDir, "Sample_TF1a_3/TF1a_3_combined_L001_fill_coordinate_sorted_q20.bam", sep = "")
              )

file.exists(bamFiles)

# for full run, use canonical chromosomes
seq_data <- sapply(bamFiles, function(x) readGAlignmentPairs(x)) # have to use readGAlignmentPairs instead of "import" function
gc()


# in next step, collapse PE alignments into "fragments", i.e. the pieces of DNA which were bound before IP
sfInit(parallel = TRUE, cpus = 6)
sfExport("seq_data")
sfExport("which")
sfLibrary(GenomicAlignments)
sfLibrary(GenomicRanges)
fragments <- sfLapply(seq_data, function(x) granges(x))
fragments <- GRangesList(fragments)
# subset for canonical chromosomes
sfExport("fragments")
fragments <- sfLapply(fragments, function(x) subsetByOverlaps(x, which))
seqlevels(fragments) <- canonicalChr
sfExport("fragments")
cov.all <- sfLapply(fragments, function(x) coverage(x))
save(cov.all, file = "/home/skurscheid/Data/Tremethick/H2A.Bbd_TF1/metagene_analysis/cov.all.rda")
save(fragments, file = "/home/skurscheid/Data/Tremethick/H2A.Bbd_TF1/metagene_analysis/fragments.rda")
sfStop()

#-----------------------------------------------------------------------------------
# BAM file processing - need to transfer into separate script
#-----------------------------------------------------------------------------------

# load
load("/home/skurscheid/Data/Tremethick/H2A.Bbd_TF1/metagene_analysis/cov.all.rda")
load("/home/skurscheid/Data/Tremethick/H2A.Bbd_TF1/metagene_analysis/fragments.rda")

# simplest and fastest way to create coverage vectors:
# 1) convert GRangesList containing all TSS into an IRangesList object
irl.tssUp1000Dn1000 <- as(grl.tssUp1000Dn1000, "IRangesList")
irl.tssUp125Dn50 <- as(grl.tssUp125Dn50, "IRangesList")
irl.intronExon <- as(grl.intronExon, "IRangesList")
# 2) use Views to resolve the RLEList into atomic vectors for each TSS region
# for all samples it should work like this - and it does
cov.tssUp1000Dn1000 <- lapply(cov.all, function(x) {v1 <- Views(x, irl.tssUp1000Dn1000); return(v1)})
cov.tssUp125Dn50 <- lapply(cov.all, function(x) {v1 <- Views(x, irl.tssUp125Dn50); return(v1)})
cov.intronExon <- lapply(cov.all, function(x) {v1 <- Views(x, irl.intronExon): return(v1)})

# Let's work with the MACS2 peaks - differential peak calling TF-1a vs TF-1

# metagene package
design <- data.frame(Samples = as.character(bamFiles),
                     align1 = c(1,1,1,0,0,0), 
                     align2 = c(0,0,0,1,1,1))
groupsFeatures <- parseFeatures(bamFiles=bamFiles, design=design, specie="human", maxDistance=1000, cores = 6)
