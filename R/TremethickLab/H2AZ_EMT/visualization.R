# load libraries
library("Gviz")
library("rtracklayer")
library("biomaRt")

# set WD
setwd("~/Data/Tremethick/MDCK_sure_select/")

# load coverage from bigwig files
h2az.wt.rep1 <- import("H2A.Z_WT_rep1/H2AZ_WT_rep1_rerun_newTrim/H2AZ-WT-rep1_S1_L001_fill_pe.bw")
h2az.wt.rep2 <- import("H2A.Z_WT_rep2/H2AZ-WT-rep2_S2_L001_fill_pe.bw")
h2az.tgfb.rep1 <- import("H2AZ_TGFb_rep1/H2AZ-TGFb-rep1_S3_L001_fill_pe.bw")
h2az.tgfb.rep2 <- import("H2AZ_TGFb_rep2/H2AZ-TGFb-rep2_S4_L001_fill_pe.bw")

