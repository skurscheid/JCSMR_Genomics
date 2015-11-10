pdf("/Volumes/gduserv/Data/Tremethick/Testes/mmus_testes/H2A_Lap1/30do_lap1_SE_NoIndex_L001_R1_001.preseq_diag.pdf")
plot(h2a.lap1.se.lc_extrap[1:90,]$TOTAL_READS / 10^6, h2a.lap1.se.lc_extrap[1:90,]$EXPECTED_DISTINCT / 10^6, type = "l", lwd = 3,
     main = "H2A.Lap1, library 30do_lap1_SE_NoIndex_R1_001",
     xlab = "Total reads (10^6)", 
     ylab = "Unique reads (10^6)")
lines(h2a.lap1.se.c_curve$total_reads / 10^6, h2a.lap1.se.c_curve$distinct_reads / 10^6, col = "red", lwd = 3)
legend("topleft", legend = c("Actual reads", "Extrapolated reads"), lwd = 3, col = c("red", "black"))
dev.off()

pdf("/Volumes/gduserv/Data/Tremethick/Testes/mmus_testes/H2A_Lap1/s_3_sequence.preseq_diag.pdf")
plot(s3.lc_extrap[1:90,]$TOTAL_READS / 10^6, s3.lc_extrap[1:90,]$EXPECTED_DISTINCT /10^6, type = "l", lwd = 3,
     main = "H2A.Lap1, library s_3_sequence.txt", 
     xlab = "Total reads (10^6)", 
     ylab = "Unique reads (10^6)")
lines(s3.c_curve$total_reads /10^6, s3.c_curve$distinct_reads / 10^6, col = "red", lwd = 3)
legend("topleft", legend = c("Actual reads", "Extrapolated reads"), lwd = 3, col = c("red", "black"))
dev.off()

tab1 <- read.table("/Volumes/gduserv/Data/Tremethick/Testes/mmus_testes/H2A_Lap1/Lap1_30_do_testis_NoIndex_L002_vs_mm10_unmasked.sorted.bam.c_curve", header = T, sep = "\t")
