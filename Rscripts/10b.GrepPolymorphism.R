#!/usr/bin/Rscript

#récupérer les positions sur les fichier filtrées selon HWE et %d'ind. génotypés:

x="un" 
pop=c("AA","BET","BRE","CEN","JAL","ODO","OIR","RIS","SAU")
  
all.files=list.files(pattern="batch_1_LF.") #, full=TRUE)
all.files2=list.files(pattern="batch_1_LP.")

 
for (k in 1:length(pop)){
	lf2=read.table(all.files[k], skip=9)
	lp2=read.table(all.files2[k], skip=9)
			
			lf2_1=lf2[,c(1:2)] #On récupère la colonne "POS" (CHROM est facultatif)
			lp2_1=lp2[,c(1:2)]
			colnames(lf2_1)=c("CHROM", "POS")
			colnames(lp2_1)=c("CHROM", "POS")
			pos_commun=merge(lf2_1, lp2_1, by=c('POS'))
			chrom=rep(x,length(pos_commun[,1]))
			position=cbind(as.character(x), pos_commun[,1])
			write.table(position, file=paste("pos_commun",pop[k], sep="."), quote=F, row.names=F, col.names=c("CHROM","POS"))
			}
 
 
 
 
 
 
 
