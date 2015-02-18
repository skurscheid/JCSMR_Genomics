# identify_LysC_Trypsin_cleavage_sites.R
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
# Purpose:
# identify LysC & Trypsin cleavage sites in peptides
# quick & dirty

# load libraries
library("gdata")
library("Biostrings")

# cardiomyocyte data
RBDpep <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/RBDmap peptides statistics.xlsx", sheet = 1, as.is = T)
RBDpep.Sequence.trimmed <- trimLRPatterns(Rpattern = "K", subject = RBDpep$RBDpep.Sequence)
RBDpep.Sequence.trimmed <- trimLRPatterns(Rpattern = "R", subject = RBDpep.Sequence.trimmed)
RBDpep.Sequence.trimmed.K <- sapply(RBDpep.Sequence.trimmed, function(x) matchPattern("K", x))
RBDpep.Sequence.trimmed.R <- sapply(RBDpep.Sequence.trimmed, function(x) matchPattern("R", x))
RBDpep.Sequence.trimmed.K.found <- unlist(lapply(RBDpep.Sequence.trimmed.K, function(x) length(x@ranges)))
RBDpep.Sequence.trimmed.R.found <- unlist(lapply(RBDpep.Sequence.trimmed.R, function(x) length(x@ranges)))
RBDpep$K_found <- as.integer(RBDpep.Sequence.trimmed.K.found)
RBDpep$R_found <- as.integer(RBDpep.Sequence.trimmed.R.found)
write.csv(RBDpep, file = "/Users/u1001407/Dropbox/REM project-Sebastian/RBDmap peptides statistics w cleavage site search.csv")

# HeLa cell data
HeLa.RBDpep <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HeLa RBDpep statistics.xlsx", sheet = 1, as.is = T)
HeLa.RBDpep <- HeLa.RBDpep[!is.na(HeLa.RBDpep$fragmentSequence),]
HeLa.RBDpep.match <- apply(HeLa.RBDpep, 1, function(x) matchPattern(x["RBDpep.Sequence"], subject = AAString(x["fragmentSequence"])))
HeLa.RBDpep.match.id <- unlist(lapply(HeLa.RBDpep.match, function(x) {x@ranges@width / x@subject@length}))
HeLa.RBDpep$Percent_ID <- HeLa.RBDpep.match.id * 100

HeLa.RBDpep.Sequence.trimmed <- trimLRPatterns(Rpattern = "K", subject = HeLa.RBDpep$RBDpep.Sequence)
HeLa.RBDpep.Sequence.trimmed <- trimLRPatterns(Rpattern = "R", subject = HeLa.RBDpep.Sequence.trimmed)
HeLa.RBDpep.Sequence.trimmed.K <- sapply(HeLa.RBDpep.Sequence.trimmed, function(x) matchPattern("K", x))
HeLa.RBDpep.Sequence.trimmed.R <- sapply(HeLa.RBDpep.Sequence.trimmed, function(x) matchPattern("R", x))
HeLa.RBDpep.Sequence.trimmed.K.found <- unlist(lapply(HeLa.RBDpep.Sequence.trimmed.K, function(x) length(x@ranges)))
HeLa.RBDpep.Sequence.trimmed.R.found <- unlist(lapply(HeLa.RBDpep.Sequence.trimmed.R, function(x) length(x@ranges)))
HeLa.RBDpep$K_found <- as.integer(HeLa.RBDpep.Sequence.trimmed.K.found)
HeLa.RBDpep$R_found <- as.integer(HeLa.RBDpep.Sequence.trimmed.R.found)
write.csv(HeLa.RBDpep, file = "/Users/u1001407/Dropbox/REM project-Sebastian/HeLa RBDpep statistics with cleavage sites.csv")
