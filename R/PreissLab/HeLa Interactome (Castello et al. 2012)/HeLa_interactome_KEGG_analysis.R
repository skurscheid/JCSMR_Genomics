#---------load libraries--------------
library("biomaRt")
library("gdata")
library("GO.db")
library("KEGGREST")
library("ggplot2")
library("gplots")
library("grid")
library("scales")

#---------custom functions------------
keggConv.batch <- function(x, max = 100, org, id.type = "ncbi-geneid") {
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
}

#---------global variables----------------------------------------
# for Fisher's Exact test
# can be "greater", "less", or "two.sided"
alternative = "greater"
p.adjust.method = "fdr"

total.keggIDs <- keggLink("hsa", "pathway")
total.keggIDs <- unique(total.keggIDs)


#---------use ENSEMBL biomaRt for annotation data-----------------
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#---------load data-------------------
kegg.brite <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/Analysis/KEGG_Brite_Hierarchy.xlsx", sheet = 1, as.is = T)
ids <- unlist(lapply(strsplit(kegg.brite$C, " "), function(x) x[1]))
rownames(kegg.brite) <- ids

#--------WCL processing----
# load data from Supplementary Information File 1 of Castello et al. 2012
wcl <- read.xls("mmc1.xls", sheet = 2)
colnames(wcl) <- c("ensembl_gene_id", "gene_symbol")

wcl.entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = wcl[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = human)
wcl.entrez_ids <- wcl.entrez_ids[-which(duplicated(wcl.entrez_ids$ensembl_gene_id)),]

wcl.entrezIDs <- unique(wcl.entrez_ids[!is.na(wcl.entrez_ids$entrezgene),]$entrezgene)
wcl.keggIDs <- keggConv.batch(wcl.entrezIDs, org = "hsa")
wcl.keggQ <- lapply(wcl.keggIDs, function(x) keggGet(x))
wcl.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(wcl.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
wcl.pathways.genes <- lapply(wcl.pathways, function(x) keggLink("genes", x))
names(wcl.pathways.genes) <- wcl.pathways
wcl.pathways.genes.entrez_ids <- unique(gsub("hsa:", "", as.character(unlist(wcl.pathways.genes))))
wcl.df <- kegg.brite[gsub("hsa", "", wcl.pathways), ]
wcl.df$ID <- rownames(wcl.df)
wcl.df$total <- rep(0, nrow(wcl.df))
wcl.df$total <- sapply(rownames(wcl.df), function(x) length(wcl.pathways.genes[[paste("hsa", x, sep = "")]]))
wcl.df$count <- rep(0, nrow(wcl.df))
wcl.df$frac <- rep(0, nrow(wcl.df))

for (i in rownames(wcl.df)) {
  kL1 <- keggLink("hsa", paste("hsa", i, sep = ""))
  wcl.df[i, ]$count <- length(which(wcl.keggIDs %in% kL1))
  wcl.df[i, ]$frac <- round(length(which(wcl.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# extract list of IDs in pathway
wcl.in_path.IDs <- lapply(rownames(wcl.df), function(x) {
  kL1 <- keggLink("hsa", paste("hsa", x, sep = ""))
  in_path <- wcl.keggIDs[which(wcl.keggIDs %in% kL1)]
})

# perform Fisher's Exact Test for each category
bkgd <- length(unique(total.keggIDs))
smpl <- length(wcl.keggIDs)
ftl <- apply(wcl.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

wcl.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
wcl.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
wcl.df$ft_fdr <- p.adjust(wcl.df$ft_pval, method = "fdr")

#--------Interactome - retrieve Entrez Gene IDs----
interactome <- read.xls("mmc1.xls", sheet = 1)
colnames(interactome)[c(1,2)] <- c("ensembl_gene_id", "gene_symbol")

# remove proteins tagged as contamination
interactome <- interactome[!interactome$mRNAbinding == "CONTAMINANT",]

interactome.entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = human)
interactome.entrez_ids <- interactome.entrez_ids[-which(duplicated(interactome.entrez_ids$ensembl_gene_id)),]
interactome.entrezIDs <- unique(interactome.entrez_ids[!is.na(interactome.entrez_ids$entrezgene),]$entrezgene)
interactome.keggIDs <- keggConv.batch(interactome.entrezIDs, org = "hsa")
interactome.keggQ <- lapply(interactome.keggIDs, function(x) keggGet(x))
interactome.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(interactome.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
interactome.pathways.genes <- lapply(interactome.pathways, function(x) keggLink("genes", x))
names(interactome.pathways.genes) <- interactome.pathways
interactome.pathways.genes.entrez_ids <- unique(gsub("hsa:", "", as.character(unlist(interactome.pathways.genes))))

interactome.df <- kegg.brite[gsub("hsa", "", interactome.pathways), ]
interactome.df$source <- rep("Interactome", nrow(interactome.df))
interactome.df$ID <- rownames(interactome.df)
# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.df), rownames(wcl.df))
interactome.df$total <- rep(0, nrow(interactome.df))
interactome.df[i1,]$total <- wcl.df[i1,]$count
interactome.df$count <- rep(0, nrow(interactome.df))
interactome.df$frac <- rep(0, nrow(interactome.df))

for (i in rownames(interactome.df)) {
  kL1 <- keggLink("hsa", paste("hsa", i, sep = ""))
  interactome.df[i, ]$count <- length(which(interactome.keggIDs %in% kL1))
  interactome.df[i, ]$frac <- round(length(which(interactome.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# extract list of IDs in pathway
interactome.in_path.IDs <- lapply(rownames(interactome.df), function(x) {
  kL1 <- keggLink("hsa", paste("hsa", x, sep = ""))
  in_path <- interactome.keggIDs[which(interactome.keggIDs %in% kL1)]
})
names(interactome.in_path.IDs) <- rownames(interactome.df)

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

# subset "no evidence"
ids <- interactome[interactome$mRNAbinding == "no evidence",]$ensembl_gene_id
interactome.ne.keggIDs <- paste("hsa:", interactome.entrez_ids[which(interactome.entrez_ids$ensembl_gene_id %in% ids),]$entrezgene, sep = "")

interactome.ne.df <- interactome.df
interactome.ne.df$source <- rep("Interactome - no evidence", nrow(interactome.ne.df))
interactome.ne.df$ID <- rownames(interactome.ne.df)
# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.ne.df), rownames(wcl.df))
interactome.ne.df$total <- 0
interactome.ne.df[i1,]$total <- wcl.df[i1,]$count
interactome.ne.df$count <- 0
interactome.ne.df$frac <- 0

for (i in rownames(interactome.ne.df)) {
  kL1 <- keggLink("hsa", paste("hsa", i, sep = ""))
  interactome.ne.df[i, ]$count <- length(which(interactome.ne.keggIDs %in% kL1))
  interactome.ne.df[i, ]$frac <- round(length(which(interactome.ne.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

ftl <- apply(interactome.ne.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.ne.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.ne.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.ne.df$ft_fdr <- p.adjust(interactome.ne.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))


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

#---------plotting------------------
df1 <- interactome.B.df[, c("A", "B", "ft_OR", "ft_fdr", "source")]
df1$ft_OR.cut <- cut(log2(df1$ft_OR), breaks = c(-2.5, 2), right = F)

dfC <- interactome.df[, c("A", "B", "C", "ft_OR", "ft_fdr", "source")]
dfC$ft_OR.cut <- cut(log2(dfC$ft_OR), breaks = c(-3.5:2), right = F)
dfC <- dfC[dfC$ft_fdr <= 0.25,]

p1 <- ggplot(data = dfC, aes(y = source, x = C)) + 
  geom_tile(aes(fill = ft_OR)) +
  theme(axis.text.y = element_text(angle = 0, size = 10), axis.title = element_blank()) +
  guides(fill = guide_legend(label.position = "bottom", direction = "horizontal")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.8, size = 12)) +
  labs(fill = "Log2 OR") 


