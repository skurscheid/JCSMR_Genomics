# LINE-1_qPCR_kallisto_analysis.R

S_H2AZ <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/S_H2AZ_vs_qPCR/abundance.txt", header = T, as.is = F)
G1_H2AZ <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/G1_H2A.Z_vs_qPCR/abundance.txt", header = T, as.is = F)
DT_H2AZ <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/DT_H2A.Z_vs_qPCR/abundance.txt", header = T, as.is = F)

S_Input <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/S_Input_vs_qPCR/abundance.txt", header = T, as.is = F)
G1_Input <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/G1_Input_vs_qPCR/abundance.txt", header = T, as.is = F)
DT_Input <- read.table("/Volumes/gduserv/Data/Tremethick/LINE_1_project/fastq/DT_Input_vs_qPCR/abundance.txt", header = T, as.is = F)

S_H2AZ.fold <- log2((S_H2AZ$tpm + 0.00001) / (S_Input$tpm + 0.00001))
G1_H2AZ.fold <- log2((G1_H2AZ$tpm + 0.00001) / (G1_Input$tpm + 0.00001))
DT_H2AZ.fold <- log2((DT_H2AZ$tpm + 0.00001) / (DT_Input$tpm + + 0.00001))

IDs <- unlist(lapply(strsplit(as.character(S_Input$target_id), "_"), function(x) x[1]))
H2AZ.fc <- data.frame(cbind(IDs, as.numeric(S_H2AZ.fold), as.numeric(G1_H2AZ.fold), as.numeric(DT_H2AZ.fold)))
colnames(H2AZ.fc)[2:4] <- c("S", "G1", "DT")

H2AZ.fc[,2] <- as.numeric(as.character(H2AZ.fc[,2]))
H2AZ.fc[,3] <- as.numeric(as.character(H2AZ.fc[,3]))
H2AZ.fc[,4] <- as.numeric(as.character(H2AZ.fc[,4]))
H2AZ.fc$mean <- apply(H2AZ.fc[,c(2:4)], 1, mean)

H2AZ.fc <- H2AZ.fc[which(H2AZ.fc$mean < 10 & H2AZ.fc$mean > -10), ]

# external function
source("~/Dropbox/Development/GeneralPurpose/R/summarySE.R")

S <- summarySE(H2AZ.fc, groupvars = "IDs", measurevar = c("S"))
G1 <- summarySE(H2AZ.fc, groupvars = "IDs", measurevar = c("G1"))
DT <- summarySE(H2AZ.fc, groupvars = "IDs", measurevar = c("DT"))
S$IDs <- paste("S", S$IDs, sep = "_")
G1$IDs <- paste("G1", G1$IDs, sep = "_")
DT$IDs <- paste("DT", DT$IDs, sep = "_")

colnames(S)[3] <- "fc"
colnames(G1)[3] <- "fc"
colnames(DT)[3] <- "fc"
J <- rbind(S, G1, DT)
J$group <- unlist(lapply(strsplit(J$IDs, "_"), function(x) x[1]))
J$group <- c(rep("S", 8), rep("G1", 8), rep("DT", 8))

p1 <- ggplot(J, aes(x = IDs, y = fc, fill = group)) 
p1 <- p1 + geom_bar(stat = "identity", position = position_dodge()) 
p1 <- p1 + facet_grid(. ~ IDs, scales = "free", space = "free")
p1 <- p1 + theme(axis.text.x = element_blank()) 
p1 <- p1 + xlab("Primers")
p1 <- p1 + ylab("[Log2 FC ChIP/Input]")
p1 <- p1 + ggtitle("H2A.Z ChIP-Seq reads,\nmean log2 fold-change ChIP/Input")
p1 <- p1 + labs(fill = "Sample")

pdf("ChIP-Seq_qPCR_quantification.pdf", paper = "a4r")
p1
dev.off()
