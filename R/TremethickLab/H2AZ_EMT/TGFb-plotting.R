# TGFb-1 plotting
grid.newpage()
vp <- viewport(width = 1, height = 1)
pushViewport(vp)

vp1 <- viewport(y = 0.75, width = 0.9, height = 0.45)
pushViewport(vp1)

gr.plot <- promoters(gr.mesenchymalMarkers, up = 1500, down = 20000)
chromosome(dT.cov.input.emt_markers.wt) <- seqnames(gr.plot)[i]
chromosome(dT.cov.input.emt_markers.tgfb) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.wt) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.tgfb) <- seqnames(gr.plot)[i]

max.y.tss <- max(max(values(dT.cov.input.emt_markers.wt)), max(values(dT.cov.input.emt_markers.tgfb)), max(values(dT.cov.h2az.emt_markers.wt)), max(values(dT.cov.h2az.emt_markers.tgfb)))
displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))

plotTracks(list(axisTrack,
                biomTrack,
                aT.ap1Sites,
                aT.nfkbSites,
                aT.TSS,
                dT.cov.input.emt_markers.wt, 
                dT.cov.h2az.emt_markers.wt,
                dT.cov.input.emt_markers.tgfb, 
                dT.cov.h2az.emt_markers.tgfb
                
),
chromosome = as(seqnames(gr.plot), "character")[i],
from = as.integer(start(gr.plot[i]), "integer"),
to = as.integer(end(gr.plot[i]), "integer"),
extend.right = 1000,
extend.left = 1000,
main = paste(gr.plot$hgnc_symbol[i], " (", width(gr.plot[i]), "bp)", sep = ""),
strand = "*",
cex.main = 0.5,
sizes = c(0.02, 0.06, 0.04, 0.04, 0.04, 0.2, 0.2, 0.2, 0.2),
scale = 0.5,
add = TRUE)

upViewport()
vp2 <- viewport(x = 0.25, y = 0.25, width = 0.5, height = 0.5)
pushViewport(vp2)
grid.rect(gp = gpar(lty="dashed"))
upViewport()
vp3 <- viewport(x = 0.75, y = 0.25, width = 0.5, height = 0.5)
pushViewport(vp3)
grid.rect(gp = gpar(lty="dashed"))
