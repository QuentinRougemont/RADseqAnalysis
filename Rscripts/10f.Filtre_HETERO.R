#!/usr/bin/Rscript



pop=c("AA","BET","BRE","CEN","JAL","ODO","OIR","RIS","SAU")
all.files=list.files(pattern="DEF.recode.vcf") #, full=TRUE)
all.files2=list.files(pattern="batch_1.sumstats")

x="un"

for (k in 1:length(pop)){
	tmp=read.table(all.files[k],skip=9)
	tmp2=read.table(all.files2[k], h=F)
	pos=tmp2[,4]+1
	het=cbind(pos,tmp2[,c(2,5,11)]) #Edit : Should compute again Het in pegas (R) and use stats from pegas because more accurate
	vcf=tmp[,c(1,2)]
	colnames(vcf)=c("CHR","pos")
	tmp3=merge(vcf,het, by="pos")
	tmp4=subset(tmp3, tmp3[,5] <= 0.50) 
	p=unique(tmp4[,1])
	posi=cbind(rep(x,length(p)),p)
	write.table(posi,file=paste(pop[k],"pos_het", sep="."), quote=F, row.names=F, col.names=c("chr","pos"))
}

#Remarque: The Het should be computed for each ecotypes separately (see my script based on Pegas in R instead)
