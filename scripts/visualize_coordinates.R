#!/user/bin/env Rscript
library(ggplot2)
#library(reshape2)

args = commandArgs(trailingOnly=TRUE)
fname = args[1]

dat=read.table(fname, skip=1)

df = data.frame(x=dat$V2, y=dat$V1, gkm_similarity = dat$V3)
dp = ggplot(df, aes(x = x, y = y, size=gkm_similarity)) + geom_point(alpha=1.0) + theme_classic() +scale_radius( range = c(0, 5), limits=c(0,0.6), breaks=seq(0,0.6,by=0.3)) +   guides(size=guide_legend()) + theme(text=element_text(size=12), legend.position="top", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) ;


ggsave(paste(fname, ".pdf", sep=""), plot=dp, width=4,height=4*1.2)

