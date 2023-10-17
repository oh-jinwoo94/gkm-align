#!/user/bin/env Rscript
library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
fname = args[1]

mat = as.matrix(read.table(fname))

# invert coordinate for inverted loci (e.g., FADS). remove the line below to visuzlie non-inverted loci. 
if((length(args) > 1) && args[2] == "invert"){
    mat <- apply(mat, 1, rev)
}

data1 <- melt(round(mat,2))

hmap=ggplot(data1, aes(x = Var1, y = Var2,fill = value))+geom_tile()+ scale_fill_continuous(low="grey", high="blue",limits=c(0,0.6), breaks=seq(0,0.6,by=0.3)) + theme_classic()+ theme(text=element_text(size=12), legend.position="top", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) + labs(fill="gkm-similarity")

ggsave(paste(fname, ".pdf", sep=""), plot=hmap, width=4,height=4*1.37)

