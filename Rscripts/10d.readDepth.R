#!/usr/bin/Rscript

pop=c("AA","BET","BRE","CEN","JAL","ODO","OIR","RIS","SAU")
all.files=list.files(pattern="coverageGLOB_FILTER.ldepth.mean") #, full=TRUE)

pdf(file="depth_median")
par(mfrow=c(3,3))

 for (k in 1:length(pop)){
	tmp=read.table(all.files[k], header=TRUE)
	median1=median(tmp[,3])
	max1=max(tmp[,3])
	min1=min(tmp[,3])
	hist(tmp[,3], breaks=30, main=paste(pop[k], "median=" , round(median1, 3), "\n", "min=", round(min1, 3), "\n", "max=", round(max1, 3),  sep=" "))
}
dev.off()

