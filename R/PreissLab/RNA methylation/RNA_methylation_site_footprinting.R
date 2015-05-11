#  RNA_methylation_site_footprinting.R
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
library(biomaRt)
library(gdata)
library(GO.db)
library(GEOquery)
library(ade4)
library(Rsamtools)
library(rtracklayer)
library(RCurl)

# AGO2 data - hg18!!!
gsm.ago2.1 <- getGEO("GSM1048187") 
gsm.ago2.2 <- getGEO("GSM1048188")

# download BED files
download.file(gsm.ago2.1@header$supplementary_file_1, "GSM1048187.bed.gz")
download.file(gsm.ago2.2@header$supplementary_file_1, "GSM1048188.bed.gz")
system("gzip -d GSM1048187.bed")
system("gzip -d GSM1048188.bed")

# AGO2 data - hg18!!!
ago2.bed1 <- import("GSM1048187.bed")
ago2.bed2 <- import("GSM1048188.bed")




