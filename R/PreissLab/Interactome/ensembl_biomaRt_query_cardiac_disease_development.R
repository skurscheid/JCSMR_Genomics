#  ensembl_biomaRt_query_cardiac_disease_development.R
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
library("biomaRt")
library("gdata")
library("GO.db")

# create objects for access to Hsap & Mmus data
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# list available filters
filters <- listFilters(mouse)

# list available attributes
attribs <- listAttributes(mouse)

pages <- attributePages(ensembl)

# Annotation obtained from the Cardiovascular Gene Ontology Annotation Initiative
# http://www.ebi.ac.uk/GOA/CVI
# accessed 2015-01-05
# ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz
# downloaded to
# /Users/skurscheid/Dropbox/REM project-Sebastian/Annotations
cvgoa <- read.table("/Users/skurscheid/Dropbox/REM project-Sebastian/Annotations/gene_association.goa_bhf-ucl", header = F, as.is = T, sep = "\t", skip = 3)

# UniProtKB records of proteins matching search "annotation:(type:disease AND (heart OR cardiac))"
# http://www.uniprot.org/uniprot/?query=annotation%3A%28type%3Adisease+AND+%28heart+OR+cardiac%29%29&sort=score
# accessed 2015-01-06
# downloaded to 
# /Users/skurscheid/Dropbox/REM project-Sebastian/Annotations
uniprot <- read.csv("/Users/skurscheid/Dropbox/REM project-Sebastian/Annotations/uniprot_search_cardiac_OR_heart_disease.csv", header = T)

# List of proteins identified to bind RNA - Yalin
interactome <- read.xls("~/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = 1)

# retrieve human Ensembl Gene IDs corresponding to UniProt proteins
ens.hsap.ids <- getBM(attributes = c("ensembl_gene_id"), values = uniprot$Entry, filters = "uniprot_swissprot", mart = human)

# Use UniProt IDs as search terms, look up Ensembl Gene IDs and retrieve IDs of mouse homologs
# UniProt list of "heart/cardiac disease" linked genes
ens.mmus.ids.uniprot <- getLDS(attributes = c("ensembl_gene_id"), filters = "uniprot_swissprot", values = uniprot$Entry, mart = human, attributesL = c("ensembl_gene_id"), martL = mouse)[,2]

# CVGOA list of genes
ens.mmus.ids.cvgoa <- getLDS(attributes = c("ensembl_gene_id"), filters = "uniprot_swissprot", values = unique(cvgoa$V2), mart = human, attributesL = c("ensembl_gene_id"), martL = mouse)[,2]

ens.mmus.cvgoa.attribs <- getBM(attributes = c("ensembl_gene_id", 
                                       "description", 
                                       "chromosome_name", 
                                       "start_position", 
                                       "end_position", 
                                       "strand", 
                                       "uniprot_swissprot", 
                                       "uniprot_genename"
                                       ),
                        filters = "ensembl_gene_id",
                        values = unique(ens.mmus.ids.cvgoa),
                        mart = mouse)

# retrieve information on homology of selected mouse genes with human genes
ens.mmus.cvgoa.attribs.homology <- getBM(attributes = c("ensembl_gene_id", 
                                       "hsapiens_homolog_ensembl_gene",
                                       "hsapiens_homolog_orthology_type",
                                       "hsapiens_homolog_orthology_confidence",
                                       "hsapiens_homolog_perc_id_r1"),
                                filters = "ensembl_gene_id",
                                values = unique(ens.mmus.ids.cvgoa),
                                mart = mouse)

# apply filter to retrieve only one2one orthologs with high confidence
ens.mmus.cvgoa.attribs.homology <- ens.mmus.cvgoa.attribs.homology[which(ens.mmus.cvgoa.attribs.homology$hsapiens_homolog_orthology_type == "ortholog_one2one" & ens.mmus.cvgoa.attribs.homology$hsapiens_homolog_orthology_confidence == 1),]
# get GO annotation for these genes
ens.mmus.cvgoa.attribs.homology.go <- getBM(attributes = c("ensembl_gene_id",
                                                           "go_id",
                                                           "goslim_goa_accession",
                                                           "goslim_goa_description"),
                                            filters = "ensembl_gene_id",
                                            values = intersect(interactome$ENSEMBL.1145, ens.mmus.cvgoa.attribs.homology$ensembl_gene_id),
                                            mart = mouse)

# only select genes with RNA binding related annotations
ens.mmus.cvgoa.attribs.homology.go.rna_bind <- ens.mmus.cvgoa.attribs.homology.go[grep("RNA", ens.mmus.cvgoa.attribs.homology.go$TERM),][grep("bind", ens.mmus.cvgoa.attribs.homology.go[grep("RNA", ens.mmus.cvgoa.attribs.homology.go$TERM),]$TERM), ]

# determine number of proteins found in whole-cell lysate which have hits in the list of mouse orthologs of the CVGOA list
WCL <- read.xls("/Users/skurscheid/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = 7)
length(intersect(ens.mmus.cvgoa.attribs.homology$ensembl_gene_id, WCL$WCL.ENSEMBL.gene.ID))
ens.mmus.cvgoa.wcl <- intersect(ens.mmus.cvgoa.attribs.homology$ensembl_gene_id, WCL$WCL.ENSEMBL.gene.ID)
length(intersect(unique(ens.mmus.ids.cvgoa), WCL$WCL.ENSEMBL.gene.ID))

# write it to TXT file 
writeLines(intersect(interactome$ENSEMBL.1145, ens.mmus.cvgoa.attribs.homology$ensembl_gene_id), "/Users/skurscheid/Dropbox/REM project-Sebastian/CVGOA_Interactome_Overlap.txt")
writeLines(intersect(ens.mmus.ids, interactome$ENSEMBL.1145), "/Users/skurscheid/Dropbox/REM project-Sebastian/UniProt_Heart_Disease_Interactome_Overlap.txt")
writeLines(unique(ens.mmus.cvgoa.attribs.homology.go.rna_bind$ensembl_gene_id), "/Users/skurscheid/Dropbox/REM project-Sebastian/CVGOA_Interactome_Overlap_Orthologs_RBP.txt")
writeLines(unique(ens.mmus.cvgoa.wcl), "/Users/skurscheid/Dropbox/REM project-Sebastian/CVGOA_WCL_Overlap.txt")


