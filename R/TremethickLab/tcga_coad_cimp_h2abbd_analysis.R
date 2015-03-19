


for (x in hm450.files) {
    print(x)
    t1 <- read.table(paste(hm450.path, x, sep = ""), header = T, as.is = T, sep = "\t")
    if (which(hm450.files == x) == 1) {
      c1 <- t1[,1]
      beta1 <- t1[,2]
      t2 <- cbind(c1,beta1)
    } else {
    t2 <- cbind(t2, t1[,2])
    }
}  
  