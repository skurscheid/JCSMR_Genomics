#  hippocampus_methylation_vs_H2A.Lap1_occupancy.R
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
#  

########################################################################
#
# Comment
# to be run on gduserve
#
########################################################################

########################################################################
#
# load libraries
#
########################################################################
library(Gviz)
library(gdata)
library(GenomicRanges)
library(rtracklayer)
library()

# working directory on gduserv
setwd("/home/skurscheid/Data/Tremethick/H2A.Lap1_paper")

# Data preparation
# Mouse hippocampus DNA methylation data obtained from:
# Weng et al 2014, BMC Neuroscience, doi:10.1186/1471-2202-15-31
# http://www.biomedcentral.com/content/supplementary/1471-2202-15-31-s3.zip
# Original XLSX files exported to CSV
# read data
bmc1 <- read.csv("BMC Neuroscience Paper/1529372002101629_add6.csv", header = T, row.names = 1, as.is = T)
bmc2 <- read.csv("BMC Neuroscience Paper/1529372002101629_add5.csv", header = T, row.names = 1, as.is = T)
bmc <- rbind(bmc1, bmc2)

# data clean up to make it comptabile with GRanges object
# replace "#" with "*" in "strand" column
bmc[which(bmc$strand == '#'),]$strand <- rep("*", length(which(bmc$strand == '#')))

# for remaining columns replace "#" with "NA"
bmc[,c(7:13)] <- apply(bmc[,c(7:13)], 2, function(x){
  x[which(x == "#")] <- rep(NA, length(which(x == "#")))
  return(x)
  }
)

gr.mm.hippo.meth <- GRanges(seqnames = bmc$chr, IRanges(start = bmc$location., end = bmc$location.), strand = bmc$strand, mcols = bmc[,c('control', 'famine', 'region', 'TSS', 'TES', 'ref_seq', 'gene', 'CGIstart', 'CGIend', 'CGIname')])
colnames(mcols(gr.mm.hippo.meth)) <- gsub("mcols.", "", colnames(mcols(gr.mm.hippo.meth)))

# import H2A.Lap1 ChIP-Seq from BigWig file
gr.mm.hippo.h2a.Lap1 <- import("Lap1_Hippo_adult7_NoIndex_L001_fill.bw")
