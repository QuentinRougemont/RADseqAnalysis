#!/usr/bin/Rscript


###################################################################################################################################
#######################Some R lines to perform part of the analysis by Rougemont et al. 2016. Mol Ecol#############################
###################################################################################################################################

library(extrafont)
#font_import()
loadfonts()

pop=c("AA","BET","BRE","CEN","JAL","ODO","OIR","RIS","SAU")
all.files=list.files(pattern=".weir.fst") #, full=TRUE)


#### Coplot Fst ########
for(i in dir(pattern="LFvsLP")){
tmp=strsplit(i, ".weir.fst",)[[1]][1] ###à changer probablemenet
tmp2=read.table(i, h=TRUE)
assign(tmp, tmp2)
 }
 

### Set negative values to zero
LFvsLP.AA[which(LFvsLP.AA[,3]<0),3]<-0
LFvsLP.BRE[which(LFvsLP.BRE[,3]<0),3]<-0
LFvsLP.BET[which(LFvsLP.BET[,3]<0),3]<-0
LFvsLP.JAL[which(LFvsLP.JAL[,3]<0),3]<-0
LFvsLP.CEN[which(LFvsLP.CEN[,3]<0),3]<-0
LFvsLP.OIR[which(LFvsLP.OIR[,3]<0),3]<-0
LFvsLP.ODO[which(LFvsLP.ODO[,3]<0),3]<-0
LFvsLP.RIS[which(LFvsLP.RIS[,3]<0),3]<-0
LFvsLP.SAU[which(LFvsLP.SAU[,3]<0),3]<-0


### Determine quantile sup for each pop
quantAA.90=quantile(LFvsLP.AA[,3], c(.90), na.rm=T) 
FstAA_sup90=subset(LFvsLP.AA, LFvsLP.AA[,3] >= quantAA.90) 
quantBET.90=quantile(LFvsLP.BET[,3], c(.90), na.rm=T) 
FstBET_sup90=subset(LFvsLP.BET, LFvsLP.BET[,3] >= quantBET.90) 
quantOIR.90=quantile(LFvsLP.OIR[,3], c(.90), na.rm=T) 
FstOIR_sup90=subset(LFvsLP.OIR, LFvsLP.OIR[,3] >= quantOIR.90) 
quantRIS.90=quantile(LFvsLP.RIS[,3], c(.90), na.rm=T) 
FstRIS_sup90=subset(LFvsLP.RIS, LFvsLP.RIS[,3] >= quantRIS.90) 

outliers_snp=merge(merge(merge(FstAA_sup90,FstBET_sup90,by="POS"),FstRIS_sup90,by="POS"),FstOIR_sup90,by="POS")
#These will constitutes the  40 SNPs used for species discrimination and hybrid identification

write.table(outliers_snp[,c(1,3,5,7,9)],"40SNPsOUTLIERS", col.names=F,row.names=F)

###################################################################################################################################

pdf(file="CoplotFST.pdf", 15,15)
par(bg="grey92")
par(mfrow=c(6,6))
for (i in 1:length(all.files))
	{
	data=read.table(all.files[i],T)
	data[which(data[,3]<0),3]<-0
	data1=merge(data[i],data[i+1], by="POS")
	plot(data1[,3],data1[,5])
	}
dev.off()
###################################################################################################################################


##Randomization pour tester si chance ou pas
nperms=1000
random_out<-"random_outliers.txt"
add.lines=F

if (!add.lines){
	write.table( cbind( "AAvsBET","AAvsOIR","AAvsRIS","BETvsOIR","BETvsRIS","OIRvsRIS"),
	file=random_out,sep=" ",
	quote=F,col.names=F,row.names=F,append=F)
}

for (i in 1:nperms)
	{

	labAA=LFvsLP.AA[sample(1:nrow(LFvsLP.AA)),2]
	labBET=LFvsLP.BET[sample(1:nrow(LFvsLP.BET)),2]
	labOIR=LFvsLP.OIR[sample(1:nrow(LFvsLP.OIR)),2]
	labRIS=LFvsLP.RIS[sample(1:nrow(LFvsLP.RIS)),2]
	AA=cbind(labAA,LFvsLP.AA[,3])
	colnames(AA)=c("POS","FST")
	BET=cbind(labBET,LFvsLP.BET[,3])
	colnames(BET)=c("POS","FST")
	OIR=cbind(labOIR,LFvsLP.OIR[,3])
	colnames(OIR)=c("POS","FST")
	RIS=cbind(labRIS,LFvsLP.RIS[,3])
	colnames(RIS)=c("POS","FST")

	quantAA=quantile(AA[,2], c(.97), na.rm=T) 
	fst.supAA=subset(AA, AA[,2] >= quantAA)
	quantBET=quantile(BET[,2], c(.97), na.rm=T) 
	fst.supBET=subset(BET, BET[,2] >= quantBET)
	quantOIR=quantile(OIR[,2], c(.97), na.rm=T)
	fst.supOIR=subset(OIR, OIR[,2] >= quantOIR) 
	quantRIS=quantile(RIS[,2], c(.97), na.rm=T) 
	fst.supRIS=subset(RIS, RIS[,2] >= quantRIS)

	AA_BET_out=merge(fst.supAA,fst.supBET,by="POS")
	AA_OIR_out=merge(fst.supAA,fst.supOIR,by="POS")
	AA_RIS_out=merge(fst.supAA,fst.supRIS,by="POS")
	BET_OIR_out=merge(fst.supBET,fst.supOIR,by="POS")
	BET_RIS_out=merge(fst.supBET,fst.supRIS,by="POS")
	OIR_RIS_out=merge(fst.supOIR,fst.supRIS,by="POS")

	write.table( cbind( dim(AA_BET_out)[1],
						dim(AA_OIR_out)[1],
						dim(AA_RIS_out)[1],
						dim(BET_OIR_out)[1],
						dim(BET_RIS_out)[1],
						dim(OIR_RIS_out)[1]),
	file=random_out,sep=" ",
	quote=F,col.names=F,row.names=F,append=T)

}

outli=read.table("random_outliers.txt",T)
w=apply(outli,2,mean)

#export table of p-value
vect=c(dim(AA)[1],dim(BB)[1],dim(CC)[1],dim(DD)[1],dim(EE)[1],dim(FF)[1])
pdf(file="null_distrib_seuil99.pdf",12,8)
par(mfrow=c(3,2))
for (i in 1:length(pouet))
	{
	hist(pouet[,i], xlim=c(0,55), xlab="distrib_outliers")
	abline(v=vect[i], lty=2, col='red')
	pval=mean(abs(pouet[,i]) >= abs(vect[i]))
	write.table(t(pval),file="probab.txt", append=T,quote=F,row.names=F, col.names=F)
	
}
dev.off()

###################################################################################################################################
#Many more tests, plots, etc should be added to this script one day.....
