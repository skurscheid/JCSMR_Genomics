
#---------load libraries---------------------------------------------
library(biomaRt)
library(gdata)
library(Biostrings)

#---------create biomaRt connections---------------------------------
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mmus.attribs <- listAttributes(mouse)


setwd("~/Data/Preiss/Interactome/Protein Abundance/")

# Comparing protein abundances between WCL and Interactom/RBDmap data
# Data being used:
# List of proteins identified in WCL:
wcl <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "WCL", as.is = T)

# List of proteins identified as mRNA interactome:
interactome <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "sheet 1" , as.is = T)

# List of proteins with RBDpep information, i.e. RBDmap list (?)
rbdpep <- read.xls("~/Dropbox/REM project-Sebastian/RBDpep analysis/HL-1 RBDmap peptide.xlsx", sheet = 1, as.is = T)


# Illumina IDs are available in biomaRt - attribute "illumina_mouseref_8_v2
wcl.IllumIDs <- getBM(c("ensembl_gene_id", "illumina_mouseref_8_v2"), filters = "ensembl_gene_id", values = wcl$WCL.ENSEMBL.gene.ID, mart = mouse)
NAs <- which(wcl.IllumIDs$illumina_mouseref_8_v2 == "")
wcl.IllumIDs <- wcl.IllumIDs[-NAs, ]
interactome.IllumIDs <- getBM(c("ensembl_gene_id", "illumina_mouseref_8_v2"), filters = "ensembl_gene_id", values = interactome$ENSEMBL.1145.most.correct, mart = mouse)
NAs <- which(interactome.IllumIDs$illumina_mouseref_8_v2 == "")
interactome.IllumIDs <- interactome.IllumIDs[-NAs, ]

IllumIDs <- unique(c(wcl.IllumIDs$illumina_mouseref_8_v2, interactome.IllumIDs$illumina_mouseref_8_v2))

# Using expression data for looking at abundance/expression levels of WCL and mRNA interactome genes, as suggested by Bernd
# https://www.evernote.com/shard/s128/nl/2147483647/9989a775-18f7-4f3c-bab4-8105feff51b9/
library(GEOquery)

# four normal HL-1 samples
gse45207 <- getGEO("GSE45207", GSEMatrix = TRUE)
gse45207.hl1_normal <- rownames(pData(gse45207[[1]])[grep("cardiomyocte untreated", pData(gse45207[[1]])[, "source_name_ch1"]),])
# get data from GEO
gse45207.hl1_normal.data <- sapply(gse45207.hl1_normal, getGEO)
# get actuall expression data from the GEOData object
gse45207.hl1_normal.data <- lapply(gse45207.hl1_normal.data, function(x) Table(x))
# collapse the list into a data.frame
gse45207.hl1_normal.data <- do.call("cbind", gse45207.hl1_normal.data)
rownames(gse45207.hl1_normal.data) <- gse45207.hl1_normal.data$GSM1099128.ID_REF
gse45207.hl1_normal.data <- gse45207.hl1_normal.data[IllumIDs,]

x <- as(unlist(gse45207.hl1_normal.data[wcl.IllumIDs$illumina_mouseref_8_v2, grep("VALUE", colnames(gse45207.hl1_normal.data))]), "numeric")
h <- density(x, kernel="gaussian")$bw # $
w <- 1 / pnorm(0, mean=x, sd=h, lower.tail=FALSE)
gse45207.hl1_normal.data.density.wcl <- density(x, bw=h, kernel="gaussian", weights=w / length(x))
gse45207.hl1_normal.data.density.wcl$y[gse45207.hl1_normal.data.density.wcl$x < 20] <- 0

x <- as(unlist(gse45207.hl1_normal.data[interactome.IllumIDs$illumina_mouseref_8_v2, grep("VALUE", colnames(gse45207.hl1_normal.data))]), "numeric")
h <- density(x, kernel="gaussian")$bw # $
w <- 1 / pnorm(0, mean=x, sd=h, lower.tail=FALSE)
gse45207.hl1_normal.data.density.interactome <- density(x, bw=h, kernel="gaussian", weights=w / length(x))
gse45207.hl1_normal.data.density.interactome$y[gse45207.hl1_normal.data.density.interactome$x < 0] <- 0

# GSE56584
gse56584 <- getGEO("GSE56584", GSEMatrix = TRUE)
gse56584.hl1_normal <- rownames(pData(gse56584[[1]])[grep("murine control serum", pData(gse56584[[1]])[, "source_name_ch1"]),])
# get data from GEO
gse56584.hl1_normal.data <- sapply(gse56584.hl1_normal, getGEO)
# get actuall expression data from the GEOData object
gse56584.hl1_normal.data <- lapply(gse56584.hl1_normal.data, function(x) Table(x))
# collapse the list into a data.frame
gse56584.hl1_normal.data <- do.call("cbind", gse56584.hl1_normal.data)
rownames(gse56584.hl1_normal.data) <- gse56584.hl1_normal.data[,1]
gse56584.hl1_normal.data <- gse56584.hl1_normal.data[IllumIDs,]
gse56584.hl1_normal.data.density.wcl <- density(as(unlist(gse56584.hl1_normal.data[wcl.IllumIDs$illumina_mouseref_8_v2, grep("VALUE", colnames(gse56584.hl1_normal.data))]), "numeric"))
gse56584.hl1_normal.data.density.interactome <- density(as(unlist(gse56584.hl1_normal.data[interactome.IllumIDs$illumina_mouseref_8_v2, grep("VALUE", colnames(gse56584.hl1_normal.data))]), "numeric"))

# plotting the two datasets
pdf("~/Dropbox/REM project-Sebastian/Revision/HL-1 gene expression.pdf", height = 8, width = 12)
par(mfrow = c(1,2))
plot(gse45207.hl1_normal.data.density.wcl,
     lwd = 5,
     xlim = c(0, 5000),
     main = "GSE45207 (Gennebaeck et al. 2013)\n normal HL-1 cells [4 samples]",
     col = "darkgrey",
     lty = 3,
     xlab = "Gene expression [cubic spline normalized, background subtracted]")
lines(gse45207.hl1_normal.data.density.interactome, lwd = 5, col = "green")
legend("topright", legend = c("WCL", "Interactome"),
       col = c("darkgrey", "green"),
       lty = c(3,1),
       lwd = 5,
       bty = "n")
plot(gse56584.hl1_normal.data.density.wcl,
     lwd = 5,
     main = "GSE56584 (Lindig et al. unpublished)\n normal HL-1 cells [4 samples]",
     col = "darkgrey",
     lty = 3,
     xlab = "Gene expression [log2 intensities, robust spline normalization]")
lines(gse56584.hl1_normal.data.density.interactome, lwd = 5, col = "green")
legend("topright", legend = c("WCL", "Interactome"),
       col = c("darkgrey", "green"),
       lty = c(3,1),
       lwd = 5,
       bty = "n")
dev.off()


