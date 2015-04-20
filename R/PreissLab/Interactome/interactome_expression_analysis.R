#  interactome_expression_analysis.R
#
#  Copyright 2015 Sebastian Kurscheid <sebastian.kurscheid@anu.edu.au>
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
library("biomaRt")
library("gdata")
library("GO.db")
library("GEOquery")
library("ade4")

mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
filters <- listFilters(mouse)
attribs <- listAttributes(mouse)

interactome <- read.xls("~/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = 1)
load("/Volumes/GDU Server SFTP/Data/Preiss/Cardiomyocyte_Interactome_expression/gpl6246.rda")
load("/Volumes/GDU Server SFTP/Data/Preiss/Cardiomyocyte_Interactome_expression/gse22292.rda")
load("/Volumes/GDU Server SFTP/Data/Preiss/Cardiomyocyte_Interactome_expression/gse49192.rda")

interactome.add_id <- getBM(attributes = c("ensembl_transcript_id", "affy_mogene_1_0_st_v1", "description", "external_gene_name"),
                            filters = "ensembl_gene_id",
                            values = unique(interactome$ENSEMBL.1145),
                            mart = mouse)

i1 <- as.character(intersect(rownames(gse22292$GSE22292_series_matrix.txt.gz@assayData$exprs), unique(interactome.add_id$affy_mogene_1_0_st_v1)))
gse22292.exprs <- gse22292$GSE22292_series_matrix.txt.gz@assayData$exprs
gse49192.exprs <- gse49192[[2]]@assayData$exprs


# selecting following samples
# Neonatal fibroblasts, re-programmed fibroblast 2wk, re-programmed fibroblast 4wk, neonatal cardiomyocyte
s1 <- rownames(gse22292[[1]]@phenoData@data[c(16,17,18,4,5,6,7,8,9,13,14,15),])

neon.fib <- rownames(gse22292[[1]]@phenoData@data[c(16,17,18),])
rep.2wk.fib <- rownames(gse22292[[1]]@phenoData@data[c(4,5,6),])
rep.4wk.fib <- rownames(gse22292[[1]]@phenoData@data[c(7,8,9),])
neon.cardiomyo <- rownames(gse22292[[1]]@phenoData@data[c(13,14,15),])

interactome.exprs <- gse22292.exprs[i1, s1]