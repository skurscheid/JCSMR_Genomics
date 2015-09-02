# TFBindingSiteSearch.R



# load libraries
library("JASPAR2014")
library("Biostrings")
library("BSgenome.Cfamiliaris.UCSC.canFam3")
library("TFBSTools")

# Using pattern matching to find putative AP-1 binding sites

ap1Motif <- DNAString("ATGAGTCAT")
genome <- BSgenome.Cfamiliaris.UCSC.canFam3
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))
m1 <- matchPattern(ap1Motif, genome$"1", max.mismatch = 1)
gr.ap1Sites <- GRanges(seqnames = "1", as(m1, "IRanges"))
subsetByOverlaps(gr.which.tss.nostrand, gr.ap1Sites)

# AP-1 is a heterodimeric protein, consisting of c-Jun/c-Fos/ATF/JDP families
# JASPAR contains matrices (all human)
# c-Jun - MA0488.1
# c-Fos - MA0476.1
# Atf1_1 - PB0004.1
# Atf1_2 - PB0108.1


# searching for NFKB sites
# NFKB1 (4 specie consensus) - MA0105.1
# NFKB1 (human) - MA0105.2
# NFKB1 (human) - MA0105.3
opts <- list()
opts[["species"]] = 9606
opts[["name"]] = "NFKB1"
opts[["all_versions"]] = TRUE
PFMatrixList <- getMatrixSet(JASPAR2014, opts = opts)
pwmList <- PWMatrixList(MA0105.2 = toPWM(PFMatrixList[["MA0105.2"]]), MA0105.3 = toPWM(PFMatrixList[["MA0105.3"]]), use.names = TRUE)
subject <- genome$"1"
sitesetList <- searchSeq(pwmList, subject, seqname = "1", min.score = "60%", strand = "*") 
gr.nfkbSites <- c(as(sitesetList[[1]], "GRanges"), as(sitesetList[[2]], "GRanges"))
     
     