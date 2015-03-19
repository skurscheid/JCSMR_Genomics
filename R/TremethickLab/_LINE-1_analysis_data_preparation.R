#  _LINE-1_analysis_data_preparation.R
#  Copyright 2015 Sebastian Kurscheid <sebastian.kurscheid@anu.edu.au>
#
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
library("GO.db")
library("GenomicFeatures")
library("GenomicAlignments")
library("rtracklayer")
library("AnnotationDbi")
library("BSgenome.Hsapiens.UCSC.hg19")
library("metagene")
library("GenomicRanges")
library("RSQLite")

# prepare biomaRt data
mm10 <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
mm10Ensembl <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

save(mm10, file = "/home/skurscheid/Data/Annotations/GRCm38_mm10/mmusculus_gene_ensembl_GRCh38.rda")
saveDb(mm10Ensembl, file = "/home/skurscheid/Data/Annotations/GRCm38_mm10/mmusculus_gene_ensembl_GRCm38_TxDB.sqlite")
# save(chromInfo, file = "/home/skurscheid/Data/Annotations/hsapiens_gene_ensembl_chromInfo_GRCh37.rda")


