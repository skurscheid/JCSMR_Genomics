library(gdata)

setwd("~/Data/Preiss/Interactome/Protein Abundance/")

# Comparing protein abundances between WCL and Interactom/RBDmap data
# Data being used:
# List of proteins identified in WCL:
wcl <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "WCL", as.is = T)

# List of proteins identified as mRNA interactome:
interactome <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "sheet 1" , as.is = T)

# List of proteins with RBDpep information, i.e. RBDmap list (?)
rbdpep <- read.xls("~/Dropbox/REM project-Sebastian/RBDpep analysis/HL-1 RBDmap peptide.xlsx", sheet = 1, as.is = T)

protein_abundance_data <- read.xls("~/Data/Preiss/Interactome/Protein Abundance/260313_cardiomyocytes 3reps and WCL.xlsx")
protein_abundance_data <- protein_abundance_data[, c("Majority.protein.IDs", "Protein.names", "mRNA.interactome", "iBAQ.Rep1", "iBAQ.Rep2", "iBAQ.Rep3", "iBAQ.WCL")]

# find RBDpep proteins in protein abundance data
g1 <- unique(unlist(sapply(rbdpep$UniProt.ID, function(x) grep(x, protein_abundance_data$Majority.protein.IDs))))

# extract kernel density of distribution of protein abundances for different groups
d.wcl <- density(protein_abundance_data[! is.na(protein_abundance_data$iBAQ.WCL),]$iBAQ.WCL)
n.wcl <- length(protein_abundance_data[! is.na(protein_abundance_data$iBAQ.WCL),]$iBAQ.WCL)

d.interactome <- density(protein_abundance_data[! is.na(protein_abundance_data$iBAQ.Rep_avg) & protein_abundance_data$mRNA.interactome == "+", "iBAQ.Rep_avg"])
n.interactome <- length(protein_abundance_data[! is.na(protein_abundance_data$iBAQ.Rep_avg) & protein_abundance_data$mRNA.interactome == "+", "iBAQ.Rep_avg"])

# discussed with Yalin, probably best to leave out.
# d.rbdpep <- density(protein_abundance_data[g1 , "iBAQ.Rep_avg"][!is.na(protein_abundance_data[g1 , "iBAQ.Rep_avg"])])
# n.rbdpep <- length(protein_abundance_data[g1 , "iBAQ.Rep_avg"][!is.na(protein_abundance_data[g1 , "iBAQ.Rep_avg"])])

# plot the distributions
pdf(paste("Protein_abundances_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".pdf", sep = ""), height = 8, width = 10)
plot(d.wcl, col = "darkgrey", xlab = "Protein abundance [iBAQ]", main = "Distributions of protein abundances",
     lty = 3, 
     lwd = 5, 
     bty = "none",
     xlim = c(2,11))
lines(d.interactome, col = "green", lwd = 5)
legend("topright", legend = c(paste("WCL [", n.wcl, "]",sep = ""), paste("Interactome [", n.interactome, "]", sep = "")),
       col = c("darkgrey", "green"),
       lty = c(3,1),
       lwd = 5,
       bty = "n")
dev.off()
