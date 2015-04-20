# ncbi_biosystems_data_preparation.R

path <- "/Volumes/MHS//workgroups/jcsmr/PreissLab/Sebastian Kurscheid//Annotations/NCBI/BioSystems/"
files <- list.files(path, pattern = ".gz")

biosystems_tables <- lapply(files, function(x) read.table(gzfile(paste(path, x, sep = "")), header = F, as.is = T, sep = "\t", fill = T, flush = T, allowEscapes = T, strip.white = T))
names(biosystems_tables) <- unlist(lapply(strsplit(files, "\\."), function(x) x[1]))

