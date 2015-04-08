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

# load libraries



library("KEGGgraph")
library("KEGG.db")
library("reactome.db")
library("biomaRt")
library("gdata")
library("GO.db")
library("KEGGREST")

source("/Users/u1001407/Dropbox/Development/JCSMR_Genomics/R/PreissLab/map_market_V2.R")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

keggConv.batch <- function(x, max = 100, org = "mmu", id.type = "ncbi-geneid") {
  if (max > 100){
    on.exit(print("Maximum number of IDs at a given time is 100"))
  } else {
    x <- paste(id.type, x, sep = ":")
    if (length(x > 100)){
      d1 <- split(x, ceiling(seq_along(x)/max))
      s1 <- sapply(d1, function(y){
         keggConv(org, y)
      })
      return(unlist(s1))
    } else {
      s1 <- sapply(x, function(y){
        keggConv(org, y)
      })
      return(unlist(s1[,1]))
    }
  }
}

# list available filters
filters <- listFilters(mouse)
# list available attributes
attribs <- listAttributes(mouse)
pages <- attributePages(mouse)

# interactome mitochondrial proteins
mitochondrial <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/Mitochondria list (184) 26March15.xlsx", sheet = 1, as.is = T)
entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = mitochondrial[,2], filters = "ensembl_gene_id", mart = mouse)
mitochondrial <- merge(mitochondrial, entrez_ids, by.x = "interactome.mito184list.1", by.y = "ensembl_gene_id", all.x = T)
colnames(mitochondrial)[4] <- "mmus_entrez_gene_id"
rm(entrez_ids)

#-----------------------------------
# total interactome
#-----------------------------------
interactome <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = 1, as.is = T)
interactome <- interactome[, c(1,2)]
colnames(interactome) <- c("gene_symbol", "ensembl_gene_id")
entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
interactome <- merge(interactome, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
interactome.entrezIDs <- unique(interactome[!is.na(interactome$entrezgene),]$entrezgene)
interactome.keggIDs <- keggConv.batch(interactome.entrezIDs)
keggQ <- lapply(interactome.keggIDs, function(x) keggGet(x))
interactome.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
interactome.pathways.genes <- lapply(interactome.pathways, function(x) keggLink("genes", x))
interactome.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(interactome.pathways.genes))))
interactome.df <- kegg.brite[gsub("mmu", "", interactome.pathways), ]
interactome.df$total <- rep(0, nrow(interactome.df))
interactome.df$total <- sapply(rownames(interactome.df), function(x) length(interactome.pathways.genes[[paste("mmu", x, sep = "")]]))
interactome.df$count <- rep(0, nrow(interactome.df))
interactome.df$frac <- rep(0, nrow(interactome.df))

for (i in rownames(interactome.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.df[i, ]$count <- length(which(interactome.keggIDs %in% kL1))
  interactome.df[i, ]$frac <- round(length(which(interactome.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

#-----------------------------------
# whole cell lysate proteome
#-----------------------------------
wcl <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "WCL", as.is = T)
colnames(wcl) <- c("gene_symbol", "ensembl_gene_id")
entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = wcl[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
wcl <- merge(wcl, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
wcl.entrezIDs <- unique(wcl[!is.na(wcl$entrezgene),]$entrezgene)
wcl.keggIDs <- keggConv.batch(wcl.entrezIDs)
keggQ <- lapply(wcl.keggIDs, function(x) keggGet(x))
wcl.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
wcl.pathways.genes <- lapply(wcl.pathways, function(x) keggLink("genes", x))
names(wcl.pathways.genes) <- wcl.pathways
wcl.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(wcl.pathways.genes))))
wcl.df <- kegg.brite[gsub("mmu", "", wcl.pathways), ]
wcl.df$total <- rep(0, nrow(wcl.df))
wcl.df$total <- sapply(rownames(wcl.df), function(x) length(wcl.pathways.genes[[paste("mmu", x, sep = "")]]))
wcl.df$count <- rep(0, nrow(wcl.df))
wcl.df$frac <- rep(0, nrow(wcl.df))

for (i in rownames(wcl.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  wcl.df[i, ]$count <- length(which(wcl.keggIDs %in% kL1))
  wcl.df[i, ]$frac <- round(length(which(wcl.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}




#------------------------------------------
mmu.eIDs <- unique(mitochondrial[!is.na(mitochondrial$mmus_entrez_gene_id),]$mmus_entrez_gene_id)
mmu.keggIDs1 <- keggConv("mmu", paste("ncbi-geneid:", mmu.eIDs[1:100], sep = ""))
mmu.keggIDs2 <- keggConv("mmu", paste("ncbi-geneid:", mmu.eIDs[101:length(eIDs)], sep = ""))
mmu.keggIDs <- c(keggIDs1, keggIDs2)

keggQ <- lapply(keggIDs, function(x) keggGet(x))
pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[2])))
pathways.genes <- lapply(pathways, function(x) keggLink("genes", x))
pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(pathways.genes))))

bsid2info <- read.csv(paste(path, "bsid2info.brief.csv", sep = ""), header = F, as.is = T)
  
linked.pathways <- lapply(pathways, function(x) {
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




