library(viridis)
library(socorro) # custom package, download like this:
                 # `devtools::install_github('ajrominger/socorro')`

# make example distributions (with high and low evenness)
xx <- seq(0, 9, length.out = 100)
pp <- ((1:length(xx)) - 0.5) / length(xx)

hiEvenDist <- dgamma(xx, 6, 2.5)
hiEvenRank <- qgamma(pp, 6, 2.5, lower.tail = FALSE)

loEvenDist <- dgamma(xx, 1, 0.4)
loEvenRank <- qgamma(pp, 1, 0.4, lower.tail = FALSE)

# colors for high and low
loCol <- viridis(3)[3]
hiCol <- viridis(3)[1]

# plotting pars
ppars <- list(mar = c(2, 2, 0, 0) + 0.5, mgp = c(1, 0, 0), cex = 1.2)


# plot density functions
pdf('~/Desktop/fig_GDE_den.pdf', width = 4, height = 4)
par(ppars)

plot(xx, hiEvenDist, ylim = range(hiEvenDist, loEvenDist),
     type = 'n', xlab = 'Genetic diversity', ylab = 'Number of OTUs',
     axes = FALSE, frame.plot = TRUE)

polygon(c(0, xx, max(xx)), c(0, hiEvenDist, 0), col = colAlpha(hiCol, 0.5),
        border = NA)
polygon(c(0, xx, max(xx)), c(0, loEvenDist, 0), col = colAlpha(loCol, 0.5),
        border = NA)

lines(xx, hiEvenDist, col = hiCol, lwd = 2)
lines(xx, loEvenDist, col = loCol, lwd = 2)

legend('topright', legend = c('High GDE', 'Low GDE'),
       pch = 22, pt.bg = colAlpha(c(hiCol, loCol), 0.5), col = c(hiCol, loCol),
       pt.lwd = 2, pt.cex = 1.75, bty = 'n')

dev.off()


# plot ranks
pdf('~/Desktop/fig_GDE_rank.pdf', width = 4, height = 4)
par(ppars)

plot(pp, hiEvenRank, ylim = range(hiEvenRank, loEvenRank), type = 'n',
     xlab = 'OTU rank', ylab = 'Genetic diversity',
     axes = FALSE, frame.plot = TRUE)

polygon(c(pp, 1, min(pp)), c(hiEvenRank, 0, 0), col = colAlpha(hiCol, 0.5),
        border = NA)
polygon(c(pp, 1, min(pp)), c(loEvenRank, 0, 0), col = colAlpha(loCol, 0.5),
        border = NA)

lines(pp, hiEvenRank, col = hiCol, lwd = 2)
lines(pp, loEvenRank, col = loCol, lwd = 2)

legend('topright', legend = c('High GDE', 'Low GDE'),
       pch = 22, pt.bg = colAlpha(c(hiCol, loCol), 0.5), col = c(hiCol, loCol),
       pt.lwd = 2, pt.cex = 1.75, bty = 'n')

dev.off()
