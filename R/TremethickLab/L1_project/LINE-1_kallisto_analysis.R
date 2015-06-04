# LINE-1_kallisto_analysis.R

S_H2AZ <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/S_H2AZ/abundance.txt", header = T, as.is = F)
G1_H2AZ <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/G1_H2AZ/abundance.txt", header = T, as.is = F)
DT_H2AZ <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/DT_H2AZ/abundance.txt", header = T, as.is = F)

S_Input <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/S_Input/abundance.txt", header = T, as.is = F)
G1_Input <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/G1_Input/abundance.txt", header = T, as.is = F)
DT_Input <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/DT_Input/abundance.txt", header = T, as.is = F)

S_H2AZ.fold <- log2((S_H2AZ$tpm + 0.00001) / (S_Input$tpm + 0.00001))
G1_H2AZ.fold <- log2((G1_H2AZ$tpm + 0.00001) / (G1_Input$tpm + 0.00001))
DT_H2AZ.fold <- log2((DT_H2AZ$tpm + 0.00001) / (DT_Input$tpm + + 0.00001))

IDs <- readLines("/Volumes//gduserv/Data/RefGenomes/mm10_GRCm38/mm10_RepeatMasker.IDs")
H2AZ.fc <- data.frame(cbind(IDs, as.numeric(S_H2AZ.fold), as.numeric(G1_H2AZ.fold), as.numeric(DT_H2AZ.fold)))

H2AZ.fc[,2] <- as.numeric(as.character(H2AZ.fc[,2]))
H2AZ.fc[,3] <- as.numeric(as.character(H2AZ.fc[,3]))
H2AZ.fc[,4] <- as.numeric(as.character(H2AZ.fc[,4]))

# L1 elements only
S <- H2AZ.fc.L1[which(H2AZ.fc.L1$V2 < 10 & H2AZ.fc.L1$V2 > 1), 2]
G1 <- H2AZ.fc.L1[which(H2AZ.fc.L1$V3 < 10 & H2AZ.fc.L1$V3 > 1),3]
DT <- H2AZ.fc.L1[which(H2AZ.fc.L1$V4 < 10 & H2AZ.fc.L1$V4 > 1), 4]
x <- c(DT, S, G1)
  
G1f <- rep("G1", length(G1))
Sf <- rep("S", length(S))
DTf <- rep("DT", length(DT))
y <- as.factor(c(DTf, Sf, G1f))
y1 <- factor(y, levels = levels(y)[c(3,2,1)])

fit1 <- lm(x ~ y1)

pdf("LINE1_quantification_boxplot.pdf", paper = "a4")

boxplot(S, G1, DT,
        axes = FALSE, frame = FALSE,
        range = 9,
        main = "H2A.Z / Input log2 fold-changes, L1 elements\n direct quantification using kallisto",
        pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),
        lwd = 2.74
        )

axis(side = 1,
     lwd = 2.75,
     at = c(1,2,3),
     labels = c(paste("S\n", "[N = ", length(S), "]\n", "Mean log2FC = ", round(mean(S), 2), sep = ""),
                paste("G1\n", "[N = ", length(G1), "]\n", "Mean log2FC = ", round(mean(G1), 2), sep = ""), 
                paste("DT **\n", "[N = ", length(DT), "]\n", "Mean log2FC = ", round(mean(DT), 2), sep = "")),
     padj = 0.8,
     cex = 1.5)

axis(side = 2, lwd = 2.75, labels = T, cex.axis = 1.5, at = seq(0,11,2), )
mtext(side = 2, "[log2 FC ChIP/Input]", line = 2.5, cex = 1.25)
   
dev.off()


# L1Md only - DT not different to S 
S <- H2AZ.fc.L1Md[which(H2AZ.fc.L1Md$V2 < 10 & H2AZ.fc.L1Md$V2 > 1), 2]
G1 <- H2AZ.fc.L1Md[which(H2AZ.fc.L1Md$V3 < 10 & H2AZ.fc.L1Md$V3 > 1),3]
DT <- H2AZ.fc.L1Md[which(H2AZ.fc.L1Md$V4 < 10 & H2AZ.fc.L1Md$V4 > 1), 4]
x <- c(DT, S, G1)

G1f <- rep("G1", length(G1))
Sf <- rep("S", length(S))
DTf <- rep("DT", length(DT))
y <- as.factor(c(DTf, Sf, G1f))
y1 <- factor(y, levels = levels(y)[c(3,2,1)])

fit1 <- lm(x ~ y1)






