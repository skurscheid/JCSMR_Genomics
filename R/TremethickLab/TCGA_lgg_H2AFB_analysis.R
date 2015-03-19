#  TCGA_lgg_H2AFB_analysis.R
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

require(vegan)
require(ade4)
require(survival)
require(coin)
require(exactRankTests)
require(matlab)
require(gtools)


i1 <- intersect(colnames(tcga.lgg.rnaseq.counts), rownames(tcga.lgg.clinical))
temp.h2afb <- cbind(tcga.lgg.rnaseq.counts[grep("H2AFB", rownames(tcga.lgg.rnaseq.counts)), i1], tcga.lgg.clinical[i1, c("status", "time", "grade")])
colnames(temp.h2afb)[1] <- "H2AFB"
temp.h2afb$time <- round(temp.h2afb$time/30, 1)


sv.grade2 <- survfit(Surv(time, status) ~ class, data = temp.h2afb[which(temp.h2afb$grade == "2"),])
sv.grade3 <- survfit(Surv(time, status) ~ class, data = temp.h2afb[which(temp.h2afb$grade == "3"),])

lsc.grade2 <- cscores(Surv(temp.h2afb[which(temp.h2afb$grade == "2") ,]$time, temp.h2afb[which(temp.h2afb$grade == "2"),]$status), int = TRUE)
lsc.grade3 <- cscores(Surv(temp.h2afb[which(temp.h2afb$grade == "3") ,]$time, temp.h2afb[which(temp.h2afb$grade == "3"),]$status), int = TRUE)

sv.grade2.PT <- perm.test(lsc.grade2 ~ class, data = temp.h2afb[which(temp.h2afb$grade == "2") ,])
sv.grade3.PT <- perm.test(lsc.grade3 ~ class, data = temp.h2afb[which(temp.h2afb$grade == "3") ,])

boxplot(temp.h2afb[tcga.lgg.samples.gradeII.complete.non_cimp,]$H2AFB, temp.h2afb[tcga.lgg.samples.gradeII.complete.cimp,]$H2AFB)
t.test(temp.h2afb[tcga.lgg.samples.gradeII.complete.non_cimp,]$H2AFB, temp.h2afb[tcga.lgg.samples.gradeII.complete.cimp,]$H2AFB)

boxplot(log2(temp.h2afb[tcga.lgg.samples.gradeII.complete.non_cimp,]$H2AFB + 1), log2(temp.h2afb[tcga.lgg.samples.gradeII.complete.cimp,]$H2AFB + 1))
t.test(log2(temp.h2afb[tcga.lgg.samples.gradeII.complete.non_cimp,]$H2AFB + 1), log2(temp.h2afb[tcga.lgg.samples.gradeII.complete.cimp,]$H2AFB + 1))

