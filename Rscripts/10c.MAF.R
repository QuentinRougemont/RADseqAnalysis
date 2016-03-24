#!/usr/bin/Rscript

#récupérer les positions sur les fichier filtrées selon HWE et %d'ind. génotypés:

x="un" 
pop=c("AA","BET","BRE","CEN","JAL","ODO","OIR","RIS","SAU")
all.files=list.files(pattern="lf.freq") #, full=TRUE)
all.files2=list.files(pattern="lp.freq")
all.files3=list.files(pattern="glob.freqglob1")

 
 for (k in 1:length(pop)){
	lf=read.table(all.files[k], skip=9)
	lp=read.table(all.files2[k], skip=9)
	glob=read.table(all.files3[k], skip=9)
	#filter based on joint maf
 				lf1=lf[lf[,5] >=0.1 & lf[,5] <= 0.9, ]
				lp1=lp[lp[,5] >=0.1 & lp[,5] <= 0.9, ]
				glob1=glob[glob[,5]>=0.05 & glob[,5] <=0.95,]
				colnames(lp1)=c("a","b","c","d","e","f","g")
				colnames(lf1)=c("a","b","c","d","e","f","g")
				colnames(glob1)=c("a","b","c","d","e","f","g")
				data_filter=rbind(lf1,lp1,glob1)
				pop_filter=unique(data_filter[,1])
				chr=rep(x,length(pop_filter))
				posi=cbind(chr,pop_filter)
	#write output:
				write.table(posi, file=paste(pop[k],"position_glob.filtered",sep="."), quote=F, row.names=F, col.names=c("chr","pos"))
}
 
 
 
 
 
