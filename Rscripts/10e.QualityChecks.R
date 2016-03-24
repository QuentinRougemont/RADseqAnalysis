#!/usr/bin/Rscript

pop=c("AA","BET","BRE","CEN","JAL","ODO","OIR","RIS","SAU")
info=c("monomorphic","biallelic","triallelic","OneSNPbyRAD","nbSNPsbyRAD")
all.files=list.files(pattern="DEF.recode.vcf") #, full=TRUE)
all.files2=list.files(pattern="batch_1.sumstats")

pdf(file="Distrib_READ_FiltreAA.pdf")
par(mfrow=c(3,3))

for (k in 1:length(pop)){
	tmp=read.table(all.files[k],skip=9)
	tmp2=read.table(all.files2[k], h=T)
	pos=tmp2[,4]+1
	sum2=cbind(pos,tmp2[,c(2,5)])
	colnames(sum2)=c("POS","LOCUS","POS_READ")
	vcf=tmp[,c(1,2)]
	colnames(vcf)=c("CHR","POS")
	t=merge(vcf,sum2)
	test1=table(t$POS_READ)
	test0=table(t$POS)
	#mean1=mean(test0)
	#max1=max(test0)
	sum1=sum(test0==1)
	sum2=sum(test0==2)
	sum3=sum(test0==3)
	y=unique(t$LOCUS) 
	z=dim(vcf)[1]/length(unique(t$LOCUS) ) #nombre de SNPs moyens par RAD
	stat=cbind(sum1,sum2,sum3,length(y),z)
	write.table(stat,"qualityCheckStat",,quote=F,append=T) #col.names=c("monomorphic","biallelic","triallelic","OneSNPbyRAD","nbSNPsbyRAD") )
	barplot(test1, main=paste(pop[k],"Distrib locus along reads",sep =" ") , xlab="position (pb)")
	abline(h=mean(test1))

}

dev.off()
