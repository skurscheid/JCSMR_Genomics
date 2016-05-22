# peptide_size_analysis.R
# Benchmarking Cardiomyocytic RBD Assignments by RBDmap Section:
# For the RBDmap data, there is little general analysis on the peptides. 
# For example, what is the size distribution and average lengths for the RBDpeps, Npeps and Xpeps? 
# If it is possible to estimate which amino acids are crosslinked or how many amino acids in the Xpeps are crosslinked on average?
# modified: 2016-03-21 - to use updated supplementary table S2 as input for calculating peptide size distributions

require(gdata)
require(Biostrings)
require(ggplot2)

setwd("~/Data/Preiss/Interactome/Revision/")

filepath <- "~/Data/Preiss/Interactome/Revision/"

xls1 <- read.xls(paste(filepath, "Table S2-Liao-R3.xlsx", sep = ""),
                 sheet = "Table S2b. RBDpep")

peptides <- xls1[, c("ENSG",
                     "trypticPeptide",
                     "proteolyticFragment",
                     "category")]
peptides <- peptides[-1,]

peptides$ENSG <- as(peptides$ENSG, "character")
peptides$trypticPeptide <- as(peptides$trypticPeptide, "character")
peptides$proteolyticFragment <- as(peptides$proteolyticFragment, "character")

aa.npep <- AAStringSet(peptides$trypticPeptide)
aa.rbdpep <- AAStringSet(peptides$proteolyticFragment)

l1 <- lapply(seq_along(1:length(aa.npep)), function(x){
  maskMotif(aa.rbdpep[[x]], aa.npep[[x]])
})

aa.xpep.width <- unlist(lapply(seq_along(1:length(aa.npep)), function(x){
  w <- width(aa.rbdpep[x]) - maskedwidth(l1[[x]])
}))
aa.xpep.width <- aa.xpep.width[!aa.xpep.width == 0]

# estimate kernel density
h <- density(aa.xpep.width, kernel="gaussian")$bw # $
dens.aa.npep <- density(width(aa.npep), bw = h)
dens.aa.rbdpep <- density(width(aa.rbdpep), bw = h)
dens.aa.xpep.width <- density(aa.xpep.width, bw = h)
w <- 1 / pnorm(0, mean=aa.xpep.width, sd=h, lower.tail=FALSE)
dens.aa.xpep.width <- density(aa.xpep.width, bw=h, kernel="gaussian", weights=w / length(aa.xpep.width))
dens.aa.xpep.width$y[dens.aa.xpep.width$x < 0] <- 0

# Plotting ----------------------------------------------------------------
pdf(paste(filepath, "Peptide_size_distribution_2016-03-21.pdf", sep = ""), paper = "a4r")
plot(dens.aa.rbdpep,
     ylim = c(0, 0.15),
     lwd = 5,
     col = "darkgrey", main = "Distribution of peptide sizes",
     xlab = "Peptide size [AA]"
     )
lines(dens.aa.npep, lwd = 5, col = "darkgreen")
lines(dens.aa.xpep.width, lwd = 5, col = "darkred")
legend("topright", legend = c("RBDpep", "Npep", "Xpep"), col = c("darkgrey", "darkgreen", "darkred"), lwd = 5)
dev.off()


# use ggplot2 for density plots ---------------------------------------------------------
aa.npep.width <- data.frame(width = width(aa.npep), pep = rep("Npep", n = length(width(aa.npep))))
aa.rbdpep.width <- data.frame(width = width(aa.rbdpep), pep = rep("RBDpep", n = length(width(aa.rbdpep))))
aa.xpep.width <- data.frame(width = aa.xpep.width, pep = rep("Xpep", n = length(aa.xpep.width)))
pep.width <- rbind(aa.npep.width, aa.rbdpep.width, aa.xpep.width)
colnames(pep.width)[2] <- "Peptide"

pdf(paste(filepath, "Peptide_size_distribution_ggplot2_2016-03-29.pdf", sep = ""), paper = "a4r")
ggplot(pep.width, aes(width, fill = Peptide)) + geom_density(alpha = 0.6, adjust = 0.95) + scale_fill_manual(values = c("darkgreen", "darkgrey", "darkred"))
dev.off()


# plotting histograms -----------------------------------------------------
xlim = c(0,170)
ylim = c(0,280)
breaks = seq(0,170,5)
lwd = 2
ylab = list("Frequency\n[count]", side = 2, cex = lwd/2, lwd = lwd, padj = -1.2)

pdf(paste(filepath, "Peptide_size_distribution_histogram_2016-03-30.pdf", sep = ""), height = 6, width = 9.6)
par(mfrow = c(3,1), mar = c(0,6,0,1))
hist(pep.width[which(pep.width$Peptide == "RBDpep"), "width"], add = F, xlim = xlim, breaks = breaks, axes = F, col = "darkgrey", main = NULL, xlab = NULL, ylab = NULL, ylim = ylim)
axis(2, lwd = lwd)
do.call(mtext, ylab)
legend("topright", legend = c("RBDpep"), col = c("darkgrey"), lwd = 5)
hist(pep.width[which(pep.width$Peptide == "Xpep"), "width"], xlim = xlim, breaks = breaks, axes = F, col = "darkred", main = NULL, xlab = NULL, ylab = NULL, ylim = ylim)
axis(2, lwd = lwd)
do.call(mtext, ylab)
legend("topright", legend = c("Xpep"), col = c("darkred"), lwd = 5)
par(mar = c(4,6,0,1))
hist(pep.width[which(pep.width$Peptide == "Npep"), "width"], add = F, xlim = xlim, breaks = breaks, axes = F, col = "darkgreen", main = NULL, xlab = NULL, ylab = NULL, ylim = ylim)
axis(2, lwd = lwd)
axis(1, at = c(0,10,20,30,40,50,100,150,170), lwd = 3)
do.call(mtext, ylab)
mtext("Peptide size [aa]", 1, lwd = 3, padj = 2.3)
legend("topright", legend = c("Npep"), col = c("darkgreen"), lwd = 5)

dev.off()


summary(pep.width[which(pep.width$Peptide == "Npep"), "width"])
