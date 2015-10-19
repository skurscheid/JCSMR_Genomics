# MSigDB data extraction
library("XML")

# downloaded msigbdb XML file from Broad Institute 2015-09-30
# parse the XML file
doc <- xmlTreeParse("~/Data/Annotations/MSigDB/msigdb_v5.0.xml", useInternalNodes = FALSE)
root <- xmlRoot(doc)
child <- xmlChildren(root)

# get STANDARD_NAME & DESCRIPTION_BRIEF attributes of all child nodes
child.standard_name <- sapply(child, function(x) xmlGetAttr(x, name = "STANDARD_NAME"))
child.description_brief <- sapply(child, function(x) xmlGetAttr(x, name = "DESCRIPTION_BRIEF"))
child.organism <- sapply(child, function(x) xmlGetAttr(x, name = "ORGANISM")) 

names(child.standard_name) <- NULL
names(child.description_brief) <- NULL
names(child.organism) <- NULL

child.df <- data.frame(cbind(as(child.standard_name, "character"), as(child.description_brief, "character"), as(child.organism, "character")))
child.df.emt <- child.df[grep("chymal", child.df[,2]),]
child.df.emt <- child.df.emt[grep("sapiens", child.df.emt[,3]),]

# manually selected
child.df[c("505", "4242", "4927", "4930", "4931", "4932", "4933", "6970", "6971", "7217", "7218", "10328"),]
MSigDB.emt_gene_id <- sapply(child[c(505, 4242, 4927:4930, 6970, 6971, 7217, 7218, 10328)], function(x) xmlGetAttr(x, name = "MEMBERS_EZID"))

# only the "Cancer Hallmark Gene Set"
MSigDB.hallmark_set_emt_gene_id <- unlist(strsplit(sapply(child[c(10328)], function(x) xmlGetAttr(x, name = "MEMBERS_EZID")), ","))
save(MSigDB.TGFb_emt_gene_id, file = "MSigDB.TGFb_emt_gene_id.rda")

# get Cfam homologs of EMT hallmark genes
MSigDB.TGFb_emt_gene_id.cfam <- getBM(attributes = c("ensembl_gene_id", "cfamiliaris_homolog_ensembl_gene"), filters = "entrezgene", values = as(MSigDB.hallmark_set_emt_gene_id, "character"), mart = human)

# geneset from microarray study of breast cancer
MSigDB.breast_emt_gene_id <- strsplit(sapply(child[c(4927:4930)], function(x) xmlGetAttr(x, name = "MEMBERS_EZID")), ",")
MSigDB.breast_emt_gene_id <- unique(unlist(MSigDB.breast_emt_gene_id))
save(MSigDB.breast_emt_gene_id, file = "MSigDB.breast_emt_gene_id.rda")

# MCF10A cells
MSigDB.mcf10a_emt_gene_id <- unique(unlist(strsplit(sapply(child[c(6970, 6971)], function(x) xmlGetAttr(x, name = "MEMBERS_EZID")), ",")))

