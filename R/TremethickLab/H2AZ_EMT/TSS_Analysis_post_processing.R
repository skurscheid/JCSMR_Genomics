# TSS analysis
# Using Eponine

tssFiles <- list.files("/Volumes/gduserv/Data/RefGenomes/Canis_familiaris/Ensembl/archive/", pattern = "*.tss.bed")
tssPath <- "/Volumes/gduserv/Data/RefGenomes/Canis_familiaris/Ensembl/archive/"
gr.tss <- unlist(GRangesList(lapply(tssFiles, function(x){
  tss <- read.table(paste(tssPath, x, sep = ""), header = F, skip = 3, as.is = T, sep = "\t")
  gr.tss <- GRanges(tss$V1, IRanges(tss$V4, tss$V5), strand = tss$V7, score = tss$V6)
})))


