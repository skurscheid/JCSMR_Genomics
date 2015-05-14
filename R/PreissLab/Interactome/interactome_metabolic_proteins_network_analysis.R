 #  interactome_metabolic_proteins_network_analysis.R
#
#  Copyright 2014 Sebastian Kurscheid <sebastian.kurscheid@anu.edu.au>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

#---------load libraries--------------
library("KEGGgraph")
library("KEGG.db")
library("reactome.db")
library("biomaRt")
library("gdata")
library("GO.db")
library("KEGGREST")
library("jsonlite")
library("ggplot2")
library("RDAVIDWebService")

source("/Users/u1001407/Dropbox/Development/GeneralPurpose/R/map_market_V2.R")

#---------custom functions------------
keggConv.batch <- function(x, max = 100, org = "mmu", id.type = "ncbi-geneid") {
  if (max > 100){
    on.exit(print("Maximum number of IDs at a given time is 100"))
  } else {
    x <- paste(id.type, x, sep = ":")
    if (length(x > 100)){
      d1 <- split(x, ceiling(seq_along(x)/max))
      s1 <- lapply(d1, function(y){
         keggConv(org, y)
      })
      return(unlist(s1))
    } else {
      d1 <- split(x, ceiling(seq_along(x)/10))
      s1 <- lapply(d1, function(y){
        keggConv(org, y)
      })
      return(unlist(s1))
    }
  }
} # TODO: edit parameters for function call 

#---------global variables----------------------------------------
# for Fisher's Exact test
# can be "greater", "less", or "two.sided"
alternative = "greater"
p.adjust.method = "fdr"

#---------use ENSEMBL biomaRt for annotation data-----------------
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# list available filters
filters <- listFilters(mouse)
# list available attributes
attribs <- listAttributes(mouse)
pages <- attributePages(mouse)
hsap.attribs <- listAttributes(human)


kegg.brite <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/Analysis/KEGG_Brite_Hierarchy.xlsx", sheet = 1, as.is = T)
ids <- unlist(lapply(strsplit(kegg.brite$C, " "), function(x) x[1]))
rownames(kegg.brite) <- ids
total.keggIDs <- keggLink("mmu", "pathway")
total.keggIDs <- unique(total.keggIDs)
length(unique(total.keggIDs))

#-------------whole cell lysate proteome------------------------------------------------------------------------------------------------
wcl <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "WCL", as.is = T)
colnames(wcl) <- c("gene_symbol", "ensembl_gene_id")
entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = wcl[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
wcl.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = wcl[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
wcl <- merge(wcl, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
wcl.entrezIDs <- unique(wcl[!is.na(wcl$entrezgene),]$entrezgene)
wcl.keggIDs <- keggConv.batch(wcl.entrezIDs)
wcl.keggQ <- lapply(wcl.keggIDs, function(x) keggGet(x))
wcl.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(wcl.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
wcl.pathways.genes <- lapply(wcl.pathways, function(x) keggLink("genes", x))
names(wcl.pathways.genes) <- wcl.pathways
wcl.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(wcl.pathways.genes))))
wcl.df <- kegg.brite[gsub("mmu", "", wcl.pathways), ]
wcl.df$ID <- rownames(wcl.df)
wcl.df$total <- rep(0, nrow(wcl.df))
wcl.df$total <- sapply(rownames(wcl.df), function(x) length(wcl.pathways.genes[[paste("mmu", x, sep = "")]]))
wcl.df$count <- rep(0, nrow(wcl.df))
wcl.df$frac <- rep(0, nrow(wcl.df))

for (i in rownames(wcl.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  wcl.df[i, ]$count <- length(which(wcl.keggIDs %in% kL1))
  wcl.df[i, ]$frac <- round(length(which(wcl.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# extract list of IDs in pathway
wcl.in_path.IDs <- lapply(rownames(wcl.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- wcl.keggIDs[which(wcl.keggIDs %in% kL1)]
})

names(wcl.in_path.IDs) <- rownames(wcl.df)

# perform Fisher's Exact Test for each category
bkgd <- length(unique(total.keggIDs))
smpl <- length(wcl.keggIDs)
ftl <- apply(wcl.df[1,], 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

wcl.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
wcl.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
wcl.df$ft_fdr <- p.adjust(wcl.df$ft_pval, method = "fdr")

#-----------total interactome----------------------------------------------------------------------------------------------------------
interactome <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "Sheet1" , as.is = T)
colnames(interactome)[c(1,2)] <- c("ensembl_gene_id", "gene_symbol")

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

interactome.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = interactome[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
interactome.mim <- getBM(attributes = c("ensembl_gene_id", "mim_morbid_accession"), values = interactome.human_homologs[,"hsapiens_homolog_ensembl_gene"], filters = "ensembl_gene_id", mart = human)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
interactome <- merge(interactome, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
rm(entrez_ids)
interactome.entrezIDs <- unique(interactome[!is.na(interactome$entrezgene),]$entrezgene)
interactome.keggIDs <- keggConv.batch(interactome.entrezIDs)
keggQ <- lapply(interactome.keggIDs, function(x) keggGet(x))
interactome.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
interactome.pathways.genes <- lapply(interactome.pathways, function(x) keggLink("genes", x))
names(interactome.pathways.genes) <- interactome.pathways
interactome.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(interactome.pathways.genes))))

# create dataframe for counting hits in pathways
interactome.df <- kegg.brite[gsub("mmu", "", interactome.pathways), ]
interactome.df$source <- rep("Interactome", nrow(interactome.df))
interactome.df$ID <- rownames(interactome.df)
# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.df), rownames(wcl.df))
interactome.df$total <- rep(0, nrow(interactome.df))
interactome.df[i1,]$total <- wcl.df[i1,]$count
interactome.df$count <- rep(0, nrow(interactome.df))
interactome.df$frac <- rep(0, nrow(interactome.df))

for (i in rownames(interactome.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.df[i, ]$count <- length(which(interactome.keggIDs %in% kL1))
  interactome.df[i, ]$frac <- round(length(which(interactome.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# extract list of IDs in pathway
interactome.in_path.IDs <- lapply(rownames(interactome.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- interactome.keggIDs[which(interactome.keggIDs %in% kL1)]
})

# perform Fisher's Exact Test for each category
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.keggIDs)

ftl <- apply(interactome.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.df$ft_fdr <- p.adjust(interactome.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#---------interactome-summarizing data at "B" level before doing Fisher's Exact test------
interactome.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.df$B))))
colnames(interactome.B.df) <- c("B", "A", "total", "count", "source")
interactome.B.df$B <- unique(interactome.df$B)
interactome.B.df$A <- sapply(unique(interactome.df$B), function(x) {A <- unique(interactome.df[which(interactome.df$B %in% x), "A"])})
interactome.B.df$source <- rep("Interactome", nrow(interactome.B.df))
interactome.B.df$total <- sapply(unique(interactome.df$B), function(x) {tot <- sum(interactome.df[which(interactome.df$B %in% x), "total"])})
interactome.B.df$count <- sapply(unique(interactome.df$B), function(x) {count <- sum(interactome.df[which(interactome.df$B %in% x), "count"])})

bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.keggIDs)

ftl <- apply(interactome.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.B.df$ft_fdr <- p.adjust(interactome.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#-----------GO RNA unrelated-------------------------------------------
interactome.go_rna_unrelated <- read.xls(xls = "/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "sheet 1", as.is = T)
colnames(interactome.go_rna_unrelated)[1:2] <- c("ensembl_gene_id", "gene_symbol")
interactome.go_rna_unrelated$GO <- as.factor(interactome.go_rna_unrelated$GO)
interactome.go_rna_unrelated$RBD <- as.factor(interactome.go_rna_unrelated$RBD)
interactome.go_rna_unrelated <- interactome.go_rna_unrelated[which(interactome.go_rna_unrelated$GO == "unrelated"),]

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome.go_rna_unrelated[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
interactome.go_rna_unrelated.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = interactome.go_rna_unrelated[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
interactome.go_rna_unrelated <- merge(interactome.go_rna_unrelated, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
interactome.go_rna_unrelated.entrezIDs <- unique(interactome.go_rna_unrelated[!is.na(interactome.go_rna_unrelated$entrezgene),]$entrezgene)
interactome.go_rna_unrelated.keggIDs <- keggConv.batch(interactome.go_rna_unrelated.entrezIDs)

# dataframe for count data
interactome.go_rna_unrelated.df <- interactome.df
interactome.go_rna_unrelated.df$source <- rep("GO_RNA_unrelated", nrow(interactome.go_rna_unrelated.df))
interactome.go_rna_unrelated.df$ID <- rownames(interactome.go_rna_unrelated.df)
interactome.go_rna_unrelated.df$total <- rep(0, nrow(interactome.go_rna_unrelated.df))

# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.go_rna_unrelated.df), rownames(wcl.df))
interactome.go_rna_unrelated.df[i1,]$total <- wcl.df[i1,]$count
interactome.go_rna_unrelated.df$count <- rep(0, nrow(interactome.go_rna_unrelated.df))

for (i in rownames(interactome.go_rna_unrelated.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.go_rna_unrelated.df[i, ]$count <- length(which(interactome.go_rna_unrelated.keggIDs %in% kL1))
}

# extract list of IDs in pathway
interactome.go_rna_unrelated.in_path.IDs <- lapply(rownames(interactome.go_rna_unrelated.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- interactome.go_rna_unrelated.keggIDs[which(interactome.go_rna_unrelated.keggIDs %in% kL1)]
})

# perform Fisher's Exact Test for each category
# Using WCL as background
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.go_rna_unrelated.keggIDs)

ftl <- apply(interactome.go_rna_unrelated.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_unrelated.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_unrelated.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.go_rna_unrelated.df$ft_fdr <- p.adjust(interactome.go_rna_unrelated.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

# summarizing data at "B" level before doing Fisher's Exact test
interactome.go_rna_unrelated.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.go_rna_unrelated.df$B))))
colnames(interactome.go_rna_unrelated.B.df) <- c("B", "A", "total", "count", "source")
interactome.go_rna_unrelated.B.df$B <- unique(interactome.go_rna_unrelated.df$B)
interactome.go_rna_unrelated.B.df$A <- sapply(unique(interactome.go_rna_unrelated.df$B), function(x) {A <- unique(interactome.go_rna_unrelated.df[which(interactome.go_rna_unrelated.df$B %in% x), "A"])})
interactome.go_rna_unrelated.B.df$source <- rep("GO_RNA_unrelated", nrow(interactome.go_rna_unrelated.B.df))
interactome.go_rna_unrelated.B.df$total <- sapply(unique(interactome.go_rna_unrelated.df$B), function(x) {tot <- sum(interactome.go_rna_unrelated.df[which(interactome.go_rna_unrelated.df$B %in% x), "total"])})
interactome.go_rna_unrelated.B.df$count <- sapply(unique(interactome.go_rna_unrelated.df$B), function(x) {count <- sum(interactome.go_rna_unrelated.df[which(interactome.go_rna_unrelated.df$B %in% x), "count"])})

# using WCL as background
ftl <- apply(interactome.go_rna_unrelated.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_unrelated.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_unrelated.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.go_rna_unrelated.B.df$ft_fdr <- p.adjust(interactome.go_rna_unrelated.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#-----------GO RNA related-------------------------------------------
interactome.go_rna_related <- read.xls(xls = "/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "sheet 1", as.is = T)
colnames(interactome.go_rna_related)[1:2] <- c("ensembl_gene_id", "gene_symbol")
interactome.go_rna_related$GO <- as.factor(interactome.go_rna_related$GO)
interactome.go_rna_related$RBD <- as.factor(interactome.go_rna_related$RBD)
interactome.go_rna_related <- interactome.go_rna_related[-which(interactome.go_rna_related$GO == "unrelated"),]

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome.go_rna_related[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
interactome.go_rna_related <- merge(interactome.go_rna_related, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
interactome.go_rna_related.entrezIDs <- unique(interactome.go_rna_related[!is.na(interactome.go_rna_related$entrezgene),]$entrezgene)
interactome.go_rna_related.keggIDs <- keggConv.batch(interactome.go_rna_related.entrezIDs)

# we are testing this subset of "interactome", therefore we include all the pathways from "interactome"
interactome.go_rna_related.df <- interactome.df
# TODO - make sure that same background is used in all tests!!!
# we are now using interactome as background to test for enrichment
i1 <- intersect(rownames(interactome.go_rna_related.df), rownames(wcl.df))
interactome.go_rna_related.df$total <- rep(0, nrow(interactome.go_rna_related.df))
interactome.go_rna_related.df[i1,]$total <- wcl.df[i1,]$count
interactome.go_rna_related.df$source <- rep("GO_RNA_related", nrow(interactome.go_rna_related.df))
interactome.go_rna_related.df$ID <- rownames(interactome.go_rna_related.df)
interactome.go_rna_related.df$count <- rep(0, nrow(interactome.go_rna_related.df))
interactome.go_rna_related.df$frac <- rep(0, nrow(interactome.go_rna_related.df))

for (i in rownames(interactome.go_rna_related.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.go_rna_related.df[i, ]$count <- length(which(interactome.go_rna_related.keggIDs %in% kL1))
}

# extract list of IDs in pathway
interactome.go_rna_related.in_path.IDs <- lapply(rownames(interactome.go_rna_related.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- interactome.go_rna_related.keggIDs[which(interactome.go_rna_related.keggIDs %in% kL1)]
})

# perform Fisher's Exact Test for each category
# Using WCL as background
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.go_rna_related.keggIDs)

ftl <- apply(interactome.go_rna_related.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_related.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_related.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
# changed p.adjust.method to "p.adjust.method" for more conservative control of p values, and set number to number of pathways in WCL
interactome.go_rna_related.df$ft_fdr <- p.adjust(interactome.go_rna_related.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

# summarizing data at "B" level before doing Fisher's Exact test
interactome.go_rna_related.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.go_rna_related.df$B))))
colnames(interactome.go_rna_related.B.df) <- c("B", "A", "total", "count", "source")
interactome.go_rna_related.B.df$B <- unique(interactome.go_rna_related.df$B)
interactome.go_rna_related.B.df$A <- sapply(unique(interactome.go_rna_related.df$B), function(x) {A <- unique(interactome.go_rna_related.df[which(interactome.go_rna_related.df$B %in% x), "A"])})
interactome.go_rna_related.B.df$source <- rep("GO_RNA_related", nrow(interactome.go_rna_related.B.df))
interactome.go_rna_related.B.df$total <- sapply(unique(interactome.go_rna_related.df$B), function(x) {tot <- sum(interactome.go_rna_related.df[which(interactome.go_rna_related.df$B %in% x), "total"])})
interactome.go_rna_related.B.df$count <- sapply(unique(interactome.go_rna_related.df$B), function(x) {count <- sum(interactome.go_rna_related.df[which(interactome.go_rna_related.df$B %in% x), "count"])})

ftl <- apply(interactome.go_rna_related.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_related.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_related.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
# changed p.adjust.method to "p.adjust.method" for more conservative control of p values, and set number to number of pathways in WCL
interactome.go_rna_related.B.df$ft_fdr <- p.adjust(interactome.go_rna_related.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#-----------plotting of KEGG enrichment analysis results---------------------------------
df1 <- rbind(interactome.B.df[, c("A", "B", "ft_OR", "ft_fdr", "source")],
             interactome.go_rna_related.B.df[, c("A", "B", "ft_OR", "ft_fdr", "source")], 
             interactome.go_rna_unrelated.B.df[, c("A", "B", "ft_OR", "ft_fdr", "source")]
             )

df1$source <- as.factor(df1$source)
df1$source <- factor(df1$source, levels = levels(df1$source)[c(3,1,2)])

df1$ft_OR.cut <- cut(log2(df1$ft_OR), breaks = c(-Inf,-4:4), right = F)

p3 <- ggplot(df1, aes(B, source)) + geom_tile(aes(fill = (df1$ft_OR.cut)))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90))
p3 <- p3 + scale_fill_gradientn(colours= c("red", "green"))
p3

# subsetting for FDR <= 0.025
c1 <- unique(as.character(df1[which(df1$ft_fdr <= 0.025),]$B))
df2 <- df1[which(df1$B %in% c1),]
df2$B <- as.factor(as.character(df2$B))
df2$ft_OR.cut <- cut(log2(df2$ft_OR), breaks = c(-Inf,-4:4), right = F)
df2$B <- factor(df2$B, levels = levels(df2$B)[df2[df2$source == "Interactome", "B"][order(df2[which(df2$source == "Interactome"),]$ft_OR.cut)]])
ggplot(data = df2, aes(x = source, y = B)) + geom_tile(aes(fill = ft_OR.cut), colour = "white") + scale_fill_brewer(palette = "PRGn") + theme(axis.text.x = element_text(angle = 90))

# subsetting Metabolism and Genetic Information Processing
df.metab <- df1[df1$A == "Metabolism",]
df.gip <- df1[df1$A == "Genetic information Processin",]
df3 <- rbind(df.metab[which(df.metab$ft_fdr <= 0.1),], df.gip[df.gip$ft_fdr <= 0.05,])
df3 <- rbind(df.metab, df.gip)
df3$B <- factor(df3$B, levels = levels(df3$B)[df2[df3$source == "Interactome", "B"][order(df2[which(df2$source == "Interactome"),]$ft_OR.cut)]])
p4 <- ggplot(df3, aes(source, B)) + geom_tile(aes(fill = df3$ft_OR.cut))
p4 <- p4 + theme(axis.text.x = element_text(angle = 90))
p4 <- p4 + scale_fill_brewer(palette = "PRGn") 
p4

#---------Plot at KEGG C level-------------------------
dfC <- rbind(interactome.df[, c("A", "B", "C", "ft_OR", "ft_fdr", "source")],
             interactome.go_rna_related.df[, c("A", "B", "C", "ft_OR", "ft_fdr", "source")], 
             interactome.go_rna_unrelated.df[, c("A", "B", "C", "ft_OR", "ft_fdr", "source")]
)

dfC$source <- as.factor(dfC$source)
dfC$source <- factor(dfC$source, levels = levels(dfC$source)[c(3,1,2)])

select1 <- unique(as.character(dfC[which(dfC$ft_fdr <= 0.1),]$C))
select2 <- interactome.go_rna_unrelated.df[which(interactome.go_rna_unrelated.df$ft_OR > 1), "C"]
select3 <- interactome.go_rna_related.df[which(interactome.go_rna_related.df$ft_OR > 1), "C"]
select4 <- unique(c(select1, select2, select3))
dfC <- dfC[which(dfC$C %in% select4),]

dfC$C <- as.factor(as.character(dfC$C))
dfC$ft_OR.cut <- cut(log2(dfC$ft_OR), breaks = c(-Inf,-4:4), right = F)
dfC$C <- factor(dfC$C, levels = levels(dfC$C)[dfC[dfC$source == "Interactome", "C"][order(dfC[which(dfC$source == "Interactome"),]$ft_OR.cut)]])
ggplot(data = dfC, aes(x = source, y = C)) + geom_tile(aes(fill = ft_OR.cut), colour = "white") + scale_fill_brewer(palette = "PRGn") + theme(axis.text.x = element_text(angle = 90))

#---------DAVID Analysis------------------------------
david <- DAVIDWebService$new(email="sebastian.kurscheid@anu.edu.au")




#-----------scratch Pad-------------------------------
mmu.eIDs <- unique(mitochondrial[!is.na(mitochondrial$mmus_entrez_gene_id),]$mmus_entrez_gene_id)
mmu.keggIDs1 <- keggConv("mmu", paste("ncbi-geneid:", mmu.eIDs[1:100], sep = ""))
mmu.keggIDs2 <- keggConv("mmu", paste("ncbi-geneid:", mmu.eIDs[101:length(eIDs)], sep = ""))
mmu.keggIDs <- c(keggIDs1, keggIDs2)

keggQ <- lapply(keggIDs, function(x) keggGet(x))
pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[2])))
pathways.genes <- lapply(pathways, function(x) keggLink("genes", x))
pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(pathways.genes))))

bsid2info <- read.csv(paste(path, "bsid2info.brief.csv", sep = ""), header = F, as.is = T)
  
interactome.linked.pathways.metab <- lapply(interactome.pathways.metabolism, function(x) {
  bsid <- data.frame(bsid2info[which(bsid2info$V3 == x), c("V1", "V3", "V4")])
  df1 <- data.frame(biosystems_tables$biosystems_biosystems_linked[which(biosystems_tables$biosystems_biosystems_linked$V1 == bsid[,"V1"]),])
  if (nrow(df1) > 0){
    df2 <- bsid2info[which(bsid2info$V1 %in% df1[, "V2"]), c("V1", "V3", "V4")]
    df1 <- merge(df1, df2, by.x = "V2", by.y = "V1", all.x = T)
    df1 <- df1[, c(2,1,4,5,3)]
  }
  l <- list(bsid, df1)
  return(l)
})

names(interactome.linked.pathways.metab) <- interactome.pathways.metabolism
interactome.metab.map <- lapply(interactome.linked.pathways.metab, function(x) {
  if (nrow(x[[2]]) > 0) {
    node1 <- rep(x[[1]][,2], nrow(x[[2]]))
    node2 <- x[[2]][,3]
    df <- data.frame(node1, node2)
    return(df)
  }
  })
names(interactome.metab.map[[1]])


linked.pathways.keggIDs <- unique(unlist(lapply(linked.pathways, function(x) x[[2]][,3])))




names(linked.pathways) <- pathways

pathways.total <- unique(c(pathways, linked.pathways.keggIDs))
mx1 <- matrix(nrow = length(pathways.total), ncol = length(mmu.keggIDs))
colnames(mx1) <- mmu.eIDs
rownames(mx1) <- pathways.total

for (i in rownames(mx1)) {
  kL1 <- gsub("mmu:", "", as.character(keggLink("mmu", i)))
  colSel <- which(colnames(mx1) %in% kL1)
  v1 <- rep(length(colSel) / length(kL1) * 100, length(colSel))
  mx1[i, colSel] <- v1
}

mx1 <- data.frame(mx1)
mx1$sum_NA <- rep(NA, nrow(mx1))
  
# only the original 104 pathways contain genes from the 187 interactome mitochondrial proteins
# therefore we create a matrix of the 104 pathways to determine their interconnection

mx1.lkd.pw <- matrix(ncol = 104, nrow = 104)
rownames(mx1.lkd.pw) <- pathways
colnames(mx1.lkd.pw) <- pathways

for (x in colnames(mx1.lkd.pw)) {
  lpa <- linked.pathways[[x]][[2]]$V3.y
  r1 <- which(rownames(mx1.lkd.pw) %in% lpa)
  if (length(r1) > 0) {
    mx1.lkd.pw[r1, x] <- rep("X", length(r1))
  }
}

for (x in rownames(mx1.lkd.pw)) {
  lpa <- linked.pathways[[x]][[2]]$V3.y
  c1 <- which(colnames(mx1.lkd.pw) %in% lpa)
  if (length(c1) > 0) {
    mx1.lkd.pw[x, c1] <- rep("X", length(c1))
  }
}



for (i in rownames(df1.sub)) {
  kL1 <- keggLink("mmu", i)
  df1.sub[i, ]$count <- length(which(mmu.keggIDs %in% kL1))
  df1.sub[i, ]$frac <- round(length(which(mmu.keggIDs %in% kL1)) / length(kL1) * 100, 1)
}

for (i in rownames(df1.sub)) {
  kL1 <- keggLink("mmu", i)
  a <- c(a, kL1)
}

# counting interactome genes in each pathway
for (i in rownames(df1)) {
  kL1 <- keggLink("mmu", i)
  df1[i, ]$count <- length(which(mmu.keggIDs %in% kL1))
  df1[i, ]$frac <- round(length(which(mmu.keggIDs %in% kL1)) / length(kL1) * 100, 1)
}

a <- character()
for (i in rownames(df1[which(df1$count >= 5),])) {
  kL1 <- keggLink("mmu", i)
  a <- c(a, kL1)
}

map.market(id = df1[which(df1$count >= 5 & !df1$class == "Global and overview maps"),]$V1,
           area = df1[which(df1$count >= 5 & !df1$class == "Global and overview maps"),]$count,
           group = df1[which(df1$count >= 5 & !df1$class == "Global and overview maps"),]$class,
           color = df1[which(df1$count >= 5 & !df1$class == "Global and overview maps"),]$count,
           main = "Captured proteins in metabolic pathways",
           lab = c(TRUE, TRUE))

map.marketV2(id = df1[which(!df1$class == "Global and overview maps"),]$V1,
           area = df1[which(!df1$class == "Global and overview maps"),]$count,
           group = df1[which(!df1$class == "Global and overview maps"),]$class,
           color = log2(df1[which(!df1$class == "Global and overview maps"),]$frac + 1),
           main = "Cellular pathway",
           lab = c(TRUE, FALSE))

# plotting
map.marketV2(id = interactome.go_rna_unrelated.df$ID,
             area = interactome.go_rna_unrelated.df$count,
             group = interactome.go_rna_unrelated.df$B,
             color = log(interactome.go_rna_unrelated.df$ft_OR),
             main = "interactome.go_rna_unrelated proteome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_interactome.go_rna_unrelated.pdf", width = 15, height = 10)
map.marketV2(id = interactome.go_rna_unrelated.df$ID,
             area = interactome.go_rna_unrelated.df$total,
             group = interactome.go_rna_unrelated.df$A,
             color = log2(interactome.go_rna_unrelated.df$ft_OR),
             main = "interactome.go_rna_unrelated proteome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

map.marketV2(id = interactome.go_rna_related.df$ID,
             area = interactome.go_rna_related.df$count,
             group = interactome.go_rna_related.df$B,
             color = log(interactome.go_rna_related.df$ft_OR),
             main = "interactome.go_rna_related proteome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_interactome.go_rna_related.pdf", width = 15, height = 10)
map.marketV2(id = interactome.go_rna_related.df$ID,
             area = interactome.go_rna_related.df$total,
             group = interactome.go_rna_related.df$A,
             color = log2(interactome.go_rna_related.df$ft_OR),
             main = "interactome.go_rna_related proteome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

#-----------HeLa/HEK293/mESC unique proteomes-----------------------------------------------------------------
unique_proteomes <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", header = T, as.is = T, sheet = "interactome unique genes")

#-----------mESC-------------------------------
# TODO: use WCL as background(?)

mesc <- data.frame(unique_proteomes[,4][-which(unique_proteomes[,4] == "")])
colnames(mesc) <- "ensembl_gene_id"

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = mesc[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
mesc.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = mesc[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

mesc <- merge(mesc, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
mesc.entrezIDs <- unique(mesc[!is.na(mesc$entrezgene),]$entrezgene)
mesc.keggIDs <- keggConv.batch(mesc.entrezIDs)
mesc.keggQ <- lapply(mesc.keggIDs, function(x) keggGet(x))
mesc.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(mesc.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
mesc.pathways.genes <- lapply(mesc.pathways, function(x) keggLink("genes", x))
names(mesc.pathways.genes) <- mesc.pathways
mesc.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(mesc.pathways.genes))))
mesc.df <- kegg.brite[gsub("mmu", "", mesc.pathways), ]
mesc.df$ID <- rownames(mesc.df)
mesc.df$total <- rep(0, nrow(mesc.df))
mesc.df$total <- sapply(rownames(mesc.df), function(x) length(mesc.pathways.genes[[paste("mmu", x, sep = "")]]))
mesc.df$count <- rep(0, nrow(mesc.df))
mesc.df$frac <- rep(0, nrow(mesc.df))

for (i in rownames(mesc.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  mesc.df[i, ]$count <- length(which(mesc.keggIDs %in% kL1))
  mesc.df[i, ]$frac <- round(length(which(mesc.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

bkgd <- length(unique(interactome.keggIDs))
smpl <- length(mesc.keggIDs)

# perform Fisher's Exact Test for each category
ftl <- apply(mesc.df[1,], 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})

mesc.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
mesc.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
mesc.df$ft_fdr <- p.adjust(mesc.df$ft_pval, method = "fdr")

pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_midlevel_mesc.pdf", width = 15, height = 10)
map.marketV2(id = mesc.df$ID,
             area = mesc.df$count,
             group = mesc.df$B,
             color = log2(mesc.df$ft_OR),
             main = "mesc specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_mesc.pdf", width = 15, height = 10)
map.marketV2(id = mesc.df$ID,
             area = mesc.df$total,
             group = mesc.df$A,
             color = log2(mesc.df$ft_OR),
             main = "mesc specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

#-----------HeLa -----------------------------
# TODO: use WCL as background(?)

hela <- data.frame(unique_proteomes[,2][-which(unique_proteomes[,2] == "")])
colnames(hela) <- "ensembl_gene_id"

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = hela[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
hela.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = hela[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

hela <- merge(hela, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
hela.entrezIDs <- unique(hela[!is.na(hela$entrezgene),]$entrezgene)
hela.keggIDs <- keggConv.batch(hela.entrezIDs)
hela.keggQ <- lapply(hela.keggIDs, function(x) keggGet(x))
hela.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(hela.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
hela.pathways.genes <- lapply(hela.pathways, function(x) keggLink("genes", x))
names(hela.pathways.genes) <- hela.pathways
hela.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(hela.pathways.genes))))
hela.df <- kegg.brite[gsub("mmu", "", hela.pathways), ]
hela.df$ID <- rownames(hela.df)
hela.df$total <- rep(0, nrow(hela.df))
hela.df$total <- sapply(rownames(hela.df), function(x) length(hela.pathways.genes[[paste("mmu", x, sep = "")]]))
hela.df$count <- rep(0, nrow(hela.df))
hela.df$frac <- rep(0, nrow(hela.df))

for (i in rownames(hela.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  hela.df[i, ]$count <- length(which(hela.keggIDs %in% kL1))
  hela.df[i, ]$frac <- round(length(which(hela.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# perform Fisher's Exact Test for each category
bkgd <- length(unique(interactome.keggIDs))
smpl <- length(hela.keggIDs)

ftl <- apply(hela.df[1,], 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})

hela.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
hela.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
hela.df$ft_fdr <- p.adjust(hela.df$ft_pval, method = "fdr")

pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_midlevel_hela.pdf", width = 15, height = 10)
map.marketV2(id = hela.df$ID,
             area = hela.df$count,
             group = hela.df$B,
             color = log2(hela.df$ft_OR),
             main = "HeLa specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_hela.pdf", width = 15, height = 10)
map.marketV2(id = hela.df$ID,
             area = hela.df$total,
             group = hela.df$A,
             color = log2(hela.df$ft_OR),
             main = "HeLa specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

#-----------HEK293-------------------------------
# TODO: use WCL as background(?)

hek293 <- data.frame(unique_proteomes[,3][-which(unique_proteomes[,3] == "")])
colnames(hek293) <- "ensembl_gene_id"

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = hek293[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
hek293.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = hek293[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

hek293 <- merge(hek293, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
hek293.entrezIDs <- unique(hek293[!is.na(hek293$entrezgene),]$entrezgene)
hek293.keggIDs <- keggConv.batch(hek293.entrezIDs)
hek293.keggQ <- lapply(hek293.keggIDs, function(x) keggGet(x))
hek293.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(hek293.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
hek293.pathways.genes <- lapply(hek293.pathways, function(x) keggLink("genes", x))
names(hek293.pathways.genes) <- hek293.pathways
hek293.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(hek293.pathways.genes))))
hek293.df <- kegg.brite[gsub("mmu", "", hek293.pathways), ]
hek293.df$ID <- rownames(hek293.df)
hek293.df$total <- rep(0, nrow(hek293.df))
hek293.df$total <- sapply(rownames(hek293.df), function(x) length(hek293.pathways.genes[[paste("mmu", x, sep = "")]]))
hek293.df$count <- rep(0, nrow(hek293.df))
hek293.df$frac <- rep(0, nrow(hek293.df))

for (i in rownames(hek293.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  hek293.df[i, ]$count <- length(which(hek293.keggIDs %in% kL1))
  hek293.df[i, ]$frac <- round(length(which(hek293.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# perform Fisher's Exact Test for each category
bkgd <- length(unique(total.keggIDs))
smpl <- length(hek293.keggIDs)

ftl <- apply(hek293.df[1,], 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})

hek293.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
hek293.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
hek293.df$ft_fdr <- p.adjust(hek293.df$ft_pval, method = "fdr")

pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_midlevel_hek293.pdf", width = 15, height = 10)
map.marketV2(id = hek293.df$ID,
             area = hek293.df$count,
             group = hek293.df$B,
             color = log2(hek293.df$ft_OR),
             main = "hek293 specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_hek293.pdf", width = 15, height = 10)
map.marketV2(id = hek293.df$ID,
             area = hek293.df$total,
             group = hek293.df$A,
             color = log2(hek293.df$ft_OR),
             main = "hek293 specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

#-----------interactome mitochondrial proteins---------
mitochondrial <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/Mitochondria list 26March15.xlsx", sheet = 1, as.is = T)
colnames(mitochondrial)[2] <- "ensembl_gene_id"
entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = mitochondrial[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
mitochondrial.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = mitochondrial[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
mitochondrial.mim <- getBM(attributes = c("ensembl_gene_id", "mim_morbid_accession"), values = mitochondrial.human_homologs[,"hsapiens_homolog_ensembl_gene"], filters = "ensembl_gene_id", mart = human)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
mitochondrial <- merge(mitochondrial, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
rm(entrez_ids)
mitochondrial.entrezIDs <- unique(mitochondrial[!is.na(mitochondrial$entrezgene),]$entrezgene)
mitochondrial.keggIDs <- keggConv.batch(mitochondrial.entrezIDs)
mitochondrial.keggQ <- lapply(mitochondrial.keggIDs, function(x) keggGet(x))
mitochondrial.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(mitochondrial.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
mitochondrial.pathways.genes <- lapply(mitochondrial.pathways, function(x) keggLink("genes", x))
names(mitochondrial.pathways.genes) <- mitochondrial.pathways
mitochondrial.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(mitochondrial.pathways.genes))))
mitochondrial.df <- kegg.brite[gsub("mmu", "", mitochondrial.pathways), ]
mitochondrial.df$ID <- rownames(mitochondrial.df)
mitochondrial.df$total <- rep(0, nrow(mitochondrial.df))
# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(mitochondrial.df), rownames(wcl.df))
mitochondrial.df[i1,]$total <- wcl.df[i1,]$count
mitochondrial.df$count <- rep(0, nrow(mitochondrial.df))
mitochondrial.df$frac <- rep(0, nrow(mitochondrial.df))

for (i in rownames(mitochondrial.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  mitochondrial.df[i, ]$count <- length(which(mitochondrial.keggIDs %in% kL1))
  mitochondrial.df[i, ]$frac <- round(length(which(mitochondrial.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# perform Fisher's Exact Test for each category
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(mitochondrial.keggIDs)

ftl <- apply(mitochondrial.df[1,], 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})

mitochondrial.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
mitochondrial.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
mitochondrial.df$ft_fdr <- p.adjust(mitochondrial.df$ft_pval, method = "fdr")

map.marketV2(id = mitochondrial.df$ID,
             area = mitochondrial.df$count,
             group = mitochondrial.df$B,
             color = log2(mitochondrial.df$ft_OR),
             main = "Mitochondrial proteins\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))

pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_mitochondrial.pdf", width = 15, height = 10)
map.marketV2(id = mitochondrial.df$ID,
             area = mitochondrial.df$total,
             group = mitochondrial.df$A,
             color = log2(mitochondrial.df$ft_OR),
             main = "Mitochondrial interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

# OMIM analysis
# first, get human homologs of mouse interactome genes
mitochondrial.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = mitochondrial[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
mitochondrial.human_homologs.mim <- getBM(attributes = c("ensembl_gene_id", "mim_morbid_accession", "mim_gene_accession"), values = mitochondrial.human_homologs[,"hsapiens_homolog_ensembl_gene"], filters = "ensembl_gene_id", mart = human)

# testing for enrichment of OMIM genes/terms in the WCL proteome
bkgd <- 14874 # total number of genes in OMIM DB
smpl <- length(unique(mitochondrial.human_homologs.mim$mim_gene_accession)) # number of genes in WCL list with OMIM gene ID
tt <- 4404 # total number of OMIM terms, i.e. bkgd has #tt 
ct <- length(unique(mitochondrial.human_homologs.mim$mim_morbid_accession)) # numer of genes in WCL list with OMIM term

matrix(c(ct, tt - ct, smpl - ct, bkgd - tt - smpl + ct), 2, 2)
fisher.test(matrix(c(ct, tt - ct, smpl - ct, bkgd - tt - smpl + ct), 2, 2), alternative = "two.sided")

# how many genes have OMIM terms?
# mouse
length(unique(mitochondrial.human_homologs[which(mitochondrial.human_homologs$hsapiens_homolog_ensembl_gene %in% (unique(mitochondrial.human_homologs.mim[which(!is.na(mitochondrial.human_homologs.mim$mim_morbid_accession)),]$ensembl_gene_id))),]$ensembl_gene_id))
# human
length(unique(mitochondrial.human_homologs.mim[which(!is.na(mitochondrial.human_homologs.mim$mim_morbid_accession)),]$ensembl_gene_id))

#-----------HL-1 unique proteome-----------------------------------------------------------------------------------------------
# TODO: use WCL as background(?)
hl1 <- data.frame(unique_proteomes[,1])
colnames(hl1) <- "ensembl_gene_id"

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = hl1[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
hl1.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = hl1[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

hl1 <- merge(hl1, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
hl1.entrezIDs <- unique(hl1[!is.na(hl1$entrezgene),]$entrezgene)
hl1.keggIDs <- keggConv.batch(hl1.entrezIDs)
hl1.keggQ <- lapply(hl1.keggIDs, function(x) keggGet(x))
hl1.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(hl1.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
hl1.pathways.genes <- lapply(hl1.pathways, function(x) keggLink("genes", x))
names(hl1.pathways.genes) <- hl1.pathways
hl1.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(hl1.pathways.genes))))
hl1.df <- kegg.brite[gsub("mmu", "", hl1.pathways), ]
hl1.df$source <- rep("HL-1", nrow(hl1.df))
# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(hl1.df), rownames(wcl.df))
hl1.df$total <- rep(0, nrow(hl1.df))
hl1.df[i1,]$total <- wcl.df[i1,]$count
hl1.df$count <- rep(0, nrow(hl1.df))
hl1.df$frac <- rep(0, nrow(hl1.df))

for (i in rownames(hl1.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  hl1.df[i, ]$count <- length(which(hl1.keggIDs %in% kL1))
  hl1.df[i, ]$frac <- round(length(which(hl1.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# perform Fisher's Exact Test for each category
bkgd <- length(unique(interactome.keggIDs))
smpl <- length(hl1.keggIDs)

ftl <- apply(hl1.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})

hl1.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
hl1.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
hl1.df$ft_fdr <- p.adjust(hl1.df$ft_pval, method = "fdr")


# summarizing data at "B" level before doing Fisher's Exact test
hl1.df.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(hl1.df$B))))
colnames(hl1.df.B.df) <- c("B", "A", "total", "count", "source")
hl1.df.B.df$B <- unique(hl1.df$B)
hl1.df.B.df$A <- sapply(unique(hl1.df$B), function(x) {A <- unique(hl1.df[which(hl1.df$B %in% x), "A"])})
hl1.df.B.df$source <- rep("HL-1", nrow(hl1.df.B.df))
hl1.df.B.df$total <- sapply(unique(hl1.df$B), function(x) {tot <- sum(hl1.df[which(hl1.df$B %in% x), "total"])})
hl1.df.B.df$count <- sapply(unique(hl1.df$B), function(x) {count <- sum(hl1.df[which(hl1.df$B %in% x), "count"])})

bkgd <- length(unique(interactome.keggIDs))
smpl <- length(hl1.keggIDs)

ftl <- apply(hl1.df.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})

hl1.df.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
hl1.df.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
hl1.df.B.df$ft_fdr <- p.adjust(hl1.df.B.df$ft_pval, method = "fdr")

# plotting
map.marketV2(id = hl1.df$ID,
             area = hl1.df$count,
             group = hl1.df$B,
             color = log(hl1.df$ft_OR),
             main = "HL-1 specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_hl1.pdf", width = 15, height = 10)
map.marketV2(id = hl1.df$ID,
             area = hl1.df$total,
             group = hl1.df$A,
             color = log2(hl1.df$ft_OR),
             main = "HL-1 specific interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()


#-----------unrecognized RBPs-------------------------
interactome.unrecognized.rbp <- read.xls(xls = "/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "Sheet1", as.is = T)
colnames(interactome.unrecognized.rbp)[1:2] <- c("ensembl_gene_id", "gene_symbol")
interactome.unrecognized.rbp$GO <- as.factor(interactome.unrecognized.rbp$GO)
interactome.unrecognized.rbp$RBD <- as.factor(interactome.unrecognized.rbp$RBD)
# tidying up category names
interactome.unrecognized.rbp[which(interactome.unrecognized.rbp$RBD == "Classical RBD"),]$RBD <- rep("classical RBD", length(which(interactome.unrecognized.rbp$RBD == "Classical RBD")))
interactome.unrecognized.rbp[which(interactome.unrecognized.rbp$RBD == "Non-Canonical RBD"),]$RBD <- rep("non classical RBD", length(which(interactome.unrecognized.rbp$RBD == "Non-Canonical RBD")))
interactome.unrecognized.rbp$RBD <- as.factor(as.character(interactome.unrecognized.rbp$RBD))
interactome.unrecognized.rbp <- interactome.unrecognized.rbp[which(interactome.unrecognized.rbp$RBD %in% c("unknown RBD")),]

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome.unrecognized.rbp[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
interactome.unrecognized.rbp.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = interactome.unrecognized.rbp[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
interactome.unrecognized.rbp <- merge(interactome.unrecognized.rbp, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
interactome.unrecognized.rbp.entrezIDs <- unique(interactome.unrecognized.rbp[!is.na(interactome.unrecognized.rbp$entrezgene),]$entrezgene)
interactome.unrecognized.rbp.keggIDs <- keggConv.batch(interactome.unrecognized.rbp.entrezIDs)

interactome.unrecognized.rbp.df <- interactome.df

# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.unrecognized.rbp.df), rownames(wcl.df))
interactome.unrecognized.rbp.df$total <- rep(0, nrow(interactome.unrecognized.rbp.df))
interactome.unrecognized.rbp.df[i1,]$total <- wcl.df[i1,]$count
interactome.unrecognized.rbp.df$source <- rep("RBD_unrecognized", nrow(interactome.unrecognized.rbp.df))
interactome.unrecognized.rbp.df$ID <- rownames(interactome.unrecognized.rbp.df)
interactome.unrecognized.rbp.df$count <- rep(0, nrow(interactome.unrecognized.rbp.df))
interactome.unrecognized.rbp.df$frac <- rep(0, nrow(interactome.unrecognized.rbp.df))

for (i in rownames(interactome.unrecognized.rbp.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.unrecognized.rbp.df[i, ]$count <- length(which(interactome.unrecognized.rbp.keggIDs %in% kL1))
}

# perform Fisher's Exact Test for each category
# Using WCL as background
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.unrecognized.rbp.keggIDs)
ftl <- apply(interactome.unrecognized.rbp.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})
interactome.unrecognized.rbp.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.unrecognized.rbp.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.unrecognized.rbp.df$ft_fdr <- p.adjust(interactome.unrecognized.rbp.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

# summarizing data at "B" level before doing Fisher's Exact test
interactome.unrecognized.rbp.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.unrecognized.rbp.df$B))))
colnames(interactome.unrecognized.rbp.B.df) <- c("B", "A", "total", "count", "source")
interactome.unrecognized.rbp.B.df$B <- unique(interactome.unrecognized.rbp.df$B)
interactome.unrecognized.rbp.B.df$A <- sapply(unique(interactome.unrecognized.rbp.df$B), function(x) {A <- unique(interactome.unrecognized.rbp.df[which(interactome.unrecognized.rbp.df$B %in% x), "A"])})
interactome.unrecognized.rbp.B.df$source <- rep("RBD_unrecognized", nrow(interactome.unrecognized.rbp.B.df))
interactome.unrecognized.rbp.B.df$total <- sapply(unique(interactome.unrecognized.rbp.df$B), function(x) {tot <- sum(interactome.unrecognized.rbp.df[which(interactome.unrecognized.rbp.df$B %in% x), "total"])})
interactome.unrecognized.rbp.B.df$count <- sapply(unique(interactome.unrecognized.rbp.df$B), function(x) {count <- sum(interactome.unrecognized.rbp.df[which(interactome.unrecognized.rbp.df$B %in% x), "count"])})

ftl <- apply(interactome.unrecognized.rbp.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})

interactome.unrecognized.rbp.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.unrecognized.rbp.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.unrecognized.rbp.B.df$ft_fdr <- p.adjust(interactome.unrecognized.rbp.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#-----------recognized RBPs-------------------------
interactome.recognized.rbp <- read.xls(xls = "/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "sheet 1", as.is = T)
colnames(interactome.recognized.rbp)[1:2] <- c("ensembl_gene_id", "gene_symbol")
interactome.recognized.rbp$GO <- as.factor(interactome.recognized.rbp$GO)
interactome.recognized.rbp$RBD <- as.factor(interactome.recognized.rbp$RBD)
# tidying up category names
interactome.recognized.rbp[which(interactome.recognized.rbp$RBD == "Classical RBD"),]$RBD <- rep("classical RBD", length(which(interactome.recognized.rbp$RBD == "Classical RBD")))
interactome.recognized.rbp[which(interactome.recognized.rbp$RBD == "Non-Canonical RBD"),]$RBD <- rep("non classical RBD", length(which(interactome.recognized.rbp$RBD == "Non-Canonical RBD")))
interactome.recognized.rbp$RBD <- as.factor(as.character(interactome.recognized.rbp$RBD))
interactome.recognized.rbp <- interactome.recognized.rbp[-which(interactome.recognized.rbp$RBD %in% c("unknown RBD")),]

entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome.recognized.rbp[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
interactome.recognized.rbp.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = interactome.recognized.rbp[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
interactome.recognized.rbp <- merge(interactome.recognized.rbp, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
interactome.recognized.rbp.entrezIDs <- unique(interactome.recognized.rbp[!is.na(interactome.recognized.rbp$entrezgene),]$entrezgene)
interactome.recognized.rbp.keggIDs <- keggConv.batch(interactome.recognized.rbp.entrezIDs)

# to test for every pathway that has been found in the interactome dataset we use that data to test against
interactome.recognized.rbp.df <- interactome.df
# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.recognized.rbp.df), rownames(wcl.df))
interactome.recognized.rbp.df$total <- rep(0, nrow(interactome.recognized.rbp.df))
interactome.recognized.rbp.df[i1,]$total <- wcl.df[i1,]$count
interactome.recognized.rbp.df$source <- rep("RBD_recognized", nrow(interactome.recognized.rbp.df))
interactome.recognized.rbp.df$ID <- rownames(interactome.recognized.rbp.df)
interactome.recognized.rbp.df$count <- rep(0, nrow(interactome.recognized.rbp.df))
interactome.recognized.rbp.df$frac <- rep(0, nrow(interactome.recognized.rbp.df))

for (i in rownames(interactome.recognized.rbp.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.recognized.rbp.df[i, ]$count <- length(which(interactome.recognized.rbp.keggIDs %in% kL1))
}

# perform Fisher's Exact Test for each category
# Using WCL as background
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.recognized.rbp.keggIDs)
ftl <- apply(interactome.recognized.rbp.df, 1, function (x) {
  ct <- as.integer(x["count"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})
interactome.recognized.rbp.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.recognized.rbp.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.recognized.rbp.df$ft_fdr <- p.adjust(interactome.recognized.rbp.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

# summarizing data at "B" level before doing Fisher's Exact test
interactome.recognized.rbp.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.recognized.rbp.df$B))))
colnames(interactome.recognized.rbp.B.df) <- c("B", "A", "total", "count", "source")
interactome.recognized.rbp.B.df$B <- unique(interactome.recognized.rbp.df$B)
interactome.recognized.rbp.B.df$A <- sapply(unique(interactome.recognized.rbp.df$B), function(x) {A <- unique(interactome.recognized.rbp.df[which(interactome.recognized.rbp.df$B %in% x), "A"])})
interactome.recognized.rbp.B.df$source <- rep("RBD_recognized", nrow(interactome.recognized.rbp.B.df))
interactome.recognized.rbp.B.df$total <- sapply(unique(interactome.recognized.rbp.df$B), function(x) {tot <- sum(interactome.recognized.rbp.df[which(interactome.recognized.rbp.df$B %in% x), "total"])})
interactome.recognized.rbp.B.df$count <- sapply(unique(interactome.recognized.rbp.df$B), function(x) {count <- sum(interactome.recognized.rbp.df[which(interactome.recognized.rbp.df$B %in% x), "count"])})

smpl <- length(interactome.recognized.rbp.keggIDs)
ftl <- apply(interactome.recognized.rbp.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = "two.sided")
})
interactome.recognized.rbp.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.recognized.rbp.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.recognized.rbp.B.df$ft_fdr <- p.adjust(interactome.recognized.rbp.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#---------total interactome rectangular tree map plotting---------------
map.marketV2(id = interactome.df$ID,
             area = interactome.df$count,
             group = interactome.df$B,
             color = log2(interactome.df$ft_OR),
             main = "Interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_interactome.pdf", width = 15, height = 10)
map.marketV2(id = interactome.df$ID,
             area = interactome.df$total,
             group = interactome.df$A,
             color = log2(interactome.df$ft_OR),
             main = "Interactome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()

#---------wcl rectangular tree map plotting---------------

map.marketV2(id = wcl.df$ID,
             area = wcl.df$count,
             group = wcl.df$B,
             color = log(wcl.df$ft_OR),
             main = "WCL proteome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))

# rectangular tree map of top level KEGG categories
pdf("/Users/u1001407/Dropbox//REM project-Sebastian/KEGG Analysis Figures/KEGG_toplevel_wcl.pdf", width = 15, height = 10)
map.marketV2(id = wcl.df$ID,
             area = wcl.df$total,
             group = wcl.df$A,
             color = log2(wcl.df$ft_OR),
             main = "WCL proteome\nKEGG pathway enrichment/depletion",
             lab = c(TRUE, FALSE))
dev.off()


