

library(gdata)
library(biomaRt)
library(GO.db)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attribs <- listAttributes(mouse)
filters <- listFilters(mouse)
attribs.hsap <- listAttributes(human)

# load IDs
cv.assoc.proteins <- read.xls("/Volumes/MHS//workgroups/jcsmr//PreissLab/Sebastian Kurscheid/Annotations//GO/cardiovascular_associated_proteins.xlsx", sheet = 1, header = T, as.is = T)
interactome <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = 1, as.is = T)
interactome <- interactome[, c(1,2)]
colnames(interactome) <- c("gene_symbol", "ensembl_gene_id")

# biomaRt attribute uniprot_swissprot
mmus.cv.assoc <- getBM(attributes = c("ensembl_gene_id", "uniprot_swissprot"), filters = "uniprot_swissprot", values = cv.assoc.proteins[which(cv.assoc.proteins$Taxon == "10090"), "ID"], mart = mouse)
hsap.cv.assoc <- getBM(attributes = c("ensembl_gene_id", "uniprot_swissprot"), filters = "uniprot_swissprot", values = cv.assoc.proteins[which(cv.assoc.proteins$Taxon == "9606"), "ID"], mart = human)
hsap.cv.assoc.mmus.homologs <- getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"), filters = "uniprot_swissprot", values = cv.assoc.proteins[which(cv.assoc.proteins$Taxon == "9606"), "ID"], mart = human)

i1 <- intersect(unique(mmus.cv.assoc$ensembl_gene_id), interactome$ensembl_gene_id)
i2 <- intersect(unique(hsap.cv.assoc.mmus.homologs$mmusculus_homolog_ensembl_gene), interactome$ensembl_gene_id)
i3 <- c(i1[which(!i1 %in% intersect(i1, i2))], i2[which(!i2 %in% intersect(i1, i2))])

interactome.go_ids <- getBM(attributes = c("ensembl_gene_id", "go_id"), filters = "ensembl_gene_id", values = interactome$ensembl_gene_id, mart = mouse)

cv.go_terms.bp <- c("GO:0007507", "GO:0048738", "GO:0008015", "GO:0050878", "GO:0001944", "GO:0042060", "GO:0006979", "GO:0016055", "GO:0006520", "GO:0050817", "GO:0006629", "GO:0006936", "GO:0048771", "GO:0051145", "GO:0007517", "GO:0042692", "GO:0048659")
cv.go_terms.cc <- c("GO:0005739", "GO:0005578")

# some plotting
xx <- as.list(GOTERM)

interactome.cv.go_terms.bp.counts <- sapply(cv.go_terms.bp, function(x) length(which(ss$go_id == x)))
df1 <- as.data.frame(interactome.cv.go_terms.bp.counts)
df1$id <- names(interactome.cv.go_terms.bp.counts)
colnames(df1) <- c("count", "id")
df1$term <- sapply(rownames(df1), function(x) xx[x][[1]]@Term)
hist1 <- ggplot(df1, aes(term, count)) + geom_histogram(stat = "identity", fill = "blue")
hist1 <- hist1 + theme(axis.text.x = element_text(angle = 90))

interactome.cv.go_terms.cc.counts <- sapply(cv.go_terms.cc, function(x) length(which(ss$go_id == x)))
df.go_cc <- as.data.frame(interactome.cv.go_terms.cc.counts)
colnames(df.go_cc)[1] <- "count"
df.go_cc$id <- rownames(df.go_cc)
df.go_cc$term <- sapply(rownames(df.go_cc), function(x) xx[x][[1]]@Term)
hist.go_cc <- ggplot(df.go_cc, aes(term, count)) + geom_histogram(stat = "identity", fill = "blue")
hist.go_cc <- hist.go_cc + theme(axis.text.x = element_text(angle = 0), plot.title = element_text("CV-associated Interactome genes\n GO CC [N = "))

pdf("/Users/u1001407/Dropbox/REM project-Sebastian/Interactome_cardiovascular_assoc_genes_GO_CC_histogram.pdf", paper = "a4")
hist.go_cc
dev.off()



