library(biomaRt)
library(XML)

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
{s <- substring(s, 2); if(strict) tolower(s) else s},
sep = "", collapse = " " )
sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


#---------use ENSEMBL biomaRt for annotation data-----------------
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# list available filters
filters <- listFilters(mouse)
# list available attributes
mmus.attribs <- listAttributes(mouse)
pages <- attributePages(mouse)
hsap.attribs <- listAttributes(human)


#---------scrape Wikipedia entry for neuronal marker names-------
html1 <- readHTMLTable("http://en.wikipedia.org/wiki/Neuronal_lineage_marker#Neural_stem_cells_markers", header = T, trim = T, as.is = T)
neuro.markers <- data.frame(html1[[2]])
neuro.markers[,1] <- as.character(neuro.markers[,1])
neuro.markers[,2] <- as.character(neuro.markers[,2])

mmus.ensembl_genes <- getBM(attributes = c("ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position", "strand"), mart = mouse)
hsap.eg1 <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"), mart = human)
hsap.eg2 <- getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"), mart = human)
hsap.ensembl_genes <- merge(hsap.eg1, hsap.eg2, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.y = T)
rm(hsap.eg1, hsap.eg2) 

#---------scan Ensembl genes for neuronal markers----------------
l2 <- lapply(neuro.markers[,2], function(x){
  markers <- unlist(strsplit(x, "; "))
  l1 <- lapply(markers, function(y){
        hgnc_symbol_rows <- c(grep(y, hsap.ensembl_genes$hgnc_symbol), grep(capwords(y, strict = T), hsap.ensembl_genes$hgnc_symbol))
        description_rows <- c(grep(y, hsap.ensembl_genes$description), grep(capwords(y, strict = T), hsap.ensembl_genes$description))
        r <- unique(c(hgnc_symbol_rows, description_rows))
        mmus_homol <- unique(hsap.ensembl_genes[r,]$mmusculus_homolog_ensembl_gene)
  })
  names(l1) <- markers
  return(l1)
})

names(l2) <- neuro.markers[,1]



