# add_Ensembl_links.R

require(gdata)

filepath <- "~/Dropbox/REM project-Sebastian/Revision/"

xls1 <- read.xls(paste(filepath, "Sebastian-list of genes.xlsx", sep = ""),
                 sheet = "list to add link")


xls1$link <- "NA"
xls1$link <- paste("http://dec2015.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", xls1$cardiomyocyte.RBPs, sep = "")
head(xls1)

write.csv(xls1, paste(filepath, "Cardiomyocyte_RBP_IDs_with_link.csv", sep = ""), quote = T)

# continue in Excel:
# Sub HyperAdd()
# 
# 'Converts each text hyperlink selected into a working hyperlink
# 
# For Each xCell In Selection
# ActiveSheet.Hyperlinks.Add Anchor:=xCell, Address:=xCell.Formula
# Next xCell
# 
# End Sub