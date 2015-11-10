# TSS analysis
# Using Eponine
chr37.tss <- read.table("/Volumes/gduserv/Data/RefGenomes/Canis_familiaris/Ensembl/archive/Canis_familiaris.CanFam3.1.dna_rm.chromosome.37.tss.bed", header = F, skip = 3, as.is = T, sep = "\t")
gr.chr37.tss <- GRanges(chr37.tss$V1, IRanges(chr37.tss$V4, chr37.tss$V5), strand = chr37.tss$V7, score = chr37.tss$V6)
