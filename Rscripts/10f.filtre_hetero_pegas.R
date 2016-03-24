##########################
# filter data based on He#
##########################
#	tmp0=as.numeric(strsplit(system("wc -l" paste=i, intern=T), " ")[[1]][1])
#	tmp=strsplit(i, "recode.vcf",)[[1]][1]


library(pegas)
pop=c("AA","BET","BRE","CEN","JAL","ODO","OIR","RIS","SAU")
 
all.files=list.files(pattern="data.LF.DEF.") #, full=TRUE)
all.files2=list.files(pattern="data.LP.DEF")

tmp0=c(15951,14527,15633,16728,17672,161115,16821,172235)

	x="un"
for (i in 1:length(pop)){
	tmp=read.vcf(file=all.files[i], from =1,to =tmp0[i])
	lf2=read.table(all.files[i], skip=9)
	tmp2=summary(as.loci(tmp))
	tmp3=sapply(tmp2, function(x) H(x$allele))
	tmp4=length(tmp3)[1]
	tmp5=rep(x, tmp4)
	tmp6=cbind(tmp5,lf2[1:tmp0[i],2],tmp3)
	pos.tmp=subset(tmp6, tmp6[,3] <=0.5)
	write.table(pos.tmp, file=paste("pos_def.PEGAS",pop[i], sep="."),, quote=F, row.names=F)
}



for(i in dir(pattern="LF.DEF")){
tmp=strsplit(i, ".recode.vcf", )[[1]][1]
tmp2=read.table(i, skip=9)
assign(tmp, tmp2)
 }
ls()


for(i in dir(pattern="LP.DEF")){
tmp=strsplit(i, ".recode.vcf", )[[1]][1]
tmp2=read.table(i, skip=9)
assign(tmp, tmp2)
 }
ls()


AA_LF=read.vcf(file= "data.LF.DEF.AA.recode.vcf", nloci = 15951, skip = 0)
BET_LF=read.vcf(file="data.LF.DEF.BET.recode.vcf", nloci = 14527, skip = 0)
BRE_LF=read.vcf(file="data.LF.DEF.BRE.recode.vcf", nloci = 15633, skip = 0)
CEN_LF=read.vcf(file="data.LF.DEF.CEN.recode.vcf", nloci = 16728, skip = 0)
JAL_LF=read.vcf(file="data.LF.DEF.JAL.recode.vcf", nloci = 17672, skip = 0)
ODO_LF=read.vcf(file="data.LF.DEF.ODO.recode.vcf", nloci = 16115, skip = 0)
OIR_LF=read.vcf(file="data.LF.DEF.OIR.recode.vcf", nloci = 16857, skip = 0)
RIS_LF=read.vcf(file="data.LF.DEF.RIS.recode.vcf", nloci = 16821, skip = 0)
SAU_LF=read.vcf(file="data.LF.DEF.SAU.recode.vcf", nloci = 172235, skip = 0)


AA_LP=read.vcf(file= "data.LP.DEF.AA.recode.vcf", nloci = 15951, skip = 0)
BET_LP=read.vcf(file="data.LP.DEF.BET.recode.vcf", nloci = 14527, skip = 0)
BRE_LP=read.vcf(file="data.LP.DEF.BRE.recode.vcf", nloci = 15633, skip = 0)
CEN_LP=read.vcf(file="data.LP.DEF.CEN.recode.vcf", nloci = 16728, skip = 0)
JAL_LP=read.vcf(file="data.LP.DEF.JAL.recode.vcf", nloci = 17672, skip = 0)
ODO_LP=read.vcf(file="data.LP.DEF.ODO.recode.vcf", nloci = 16115, skip = 0)
OIR_LP=read.vcf(file="data.LP.DEF.OIR.recode.vcf", nloci = 16857, skip = 0)
RIS_LP=read.vcf(file="data.LP.DEF.RIS.recode.vcf", nloci = 16821, skip = 0)
SAU_LP=read.vcf(file="data.LP.DEF.SAU.recode.vcf", nloci = 172235, skip = 0)



S <- summary(as.loci(AA_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
posi=cbind(chr,data.LF.DEF.AA[,2], H1)
colnames(posi)=c("CHR","POSI", "HET")
posi2=subset(posi, posi[,3] <= 0.5)
dim(posi2) 
write.table(posi2[,c(1:2)], "pos_het_lf.AA", quote=F, row.names=F)

S <- summary(as.loci(BET_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.BET[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.BET", quote=F, row.names=F)


S <- summary(as.loci(BRE_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.BRE[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.BRE", quote=F, row.names=F)


S <- summary(as.loci(JAL_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.JAL[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.JAL", quote=F, row.names=F)

S <- summary(as.loci(ODO_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.ODO[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.ODO", quote=F, row.names=F)

S <- summary(as.loci(OIR_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.OIR[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.OIR", quote=F, row.names=F)

S <- summary(as.loci(RIS_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.RIS[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.RIS", quote=F, row.names=F)

S <- summary(as.loci(SAU_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.SAU[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.SAU", quote=F, row.names=F)



######## LP ###################
S <- summary(as.loci(AA_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.AA[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.AA", quote=F, row.names=F)

S <- summary(as.loci(BET_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.BET[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.BET", quote=F, row.names=F)


S <- summary(as.loci(BRE_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.BRE[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.BRE", quote=F, row.names=F)

S <- summary(as.loci(JAL_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.JAL[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.JAL", quote=F, row.names=F)

S <- summary(as.loci(ODO_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.ODO[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.ODO", quote=F, row.names=F)

S <- summary(as.loci(OIR_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.OIR[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.OIR", quote=F, row.names=F)

S <- summary(as.loci(RIS_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.RIS[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.RIS", quote=F, row.names=F)

S <- summary(as.loci(SAU_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.SAU[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.SAU", quote=F, row.names=F)

###CENS
S <- summary(as.loci(CEN_LP))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LP.DEF.CEN[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lp.CEN", quote=F, row.names=F)


S <- summary(as.loci(CEN_LF))
## compute H for all loci:
H1=sapply(S, function(x) H(x$allele))
z=length(S)[1]
x="un"
chr=rep(x,z)
HET_O=cbind(chr,data.LF.DEF.CEN[,2], H1)
colnames(HET_O)=c("CHR","POS", "HET_O")
HET_O2=subset(HET_O, HET_O[,3] <= 0.5)
dim(HET_O2) 
write.table(HET_O2[,c(1:2)], "pos_het_lf.CEN", quote=F, row.names=F)

rm(list=ls())
for (i in dir (pattern='pos_het_lf')){
tmp=strsplit(i, "pos_het_",)[[1]][2] ###à changer probablemenet
tmp2=read.table(i, h=T)
assign(tmp, tmp2)
}

for (i in dir (pattern='pos_het_lp')){
tmp=strsplit(i, "pos_het_",)[[1]][2] ###à changer probablemenet
tmp2=read.table(i, h=T)
assign(tmp, tmp2)
}

aa=union(lf.AA[,2],lp.AA[,2])
bet=union(lf.BET[,2],lp.BET[,2])
bre=union(lf.BRE[,2],lp.BRE[,2])
cen=union(lf.CEN[,2],lp.CEN[,2])
jal=union(lf.JAL[,2],lp.JAL[,2])
odo=union(lf.ODO[,2],lp.ODO[,2])
oir=union(lf.OIR[,2],lp.OIR[,2])
ris=union(lf.RIS[,2],lp.RIS[,2])
sau=union(lf.SAU[,2],lp.SAU[,2])

x="un"  #On creer un vecteur de CHROM de même longueur que le nb de positiion
aa1=length(aa)
bet1=length(bet)
bre1=length(bre)
cen1=length(cen)
jal1=length(jal)
odo1=length(odo)
oir1=length(oir)
ris1=length(ris)
sau1=length(sau)

aa2 =rep(x,aa1 )
bet2=rep(x,bet1)
bre2=rep(x,bre1)
cen2=rep(x,cen1)
jal2=rep(x,jal1)
odo2=rep(x,odo1)
oir2=rep(x,oir1)
ris2=rep(x,ris1)
sau2=rep(x,sau1)

aa2 =cbind(aa2 , aa )
bet2=cbind(bet2, bet)
bre2=cbind(bre2, bre)
cen2=cbind(cen2, cen)
jal2=cbind(jal2, jal)
odo2=cbind(odo2, odo)
oir2=cbind(oir2, oir)
ris2=cbind(ris2, ris)
sau2=cbind(sau2, sau)


write.table(aa2 , "pos_communHET.AA ", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(bet2, "pos_communHET.BET", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(bre2, "pos_communHET.BRE", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(cen2, "pos_communHET.CEN", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(jal2, "pos_communHET.JAL", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(odo2, "pos_communHET.ODO", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(oir2, "pos_communHET.OIR", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(ris2, "pos_communHET.RIS", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools
write.table(sau2, "pos_communHET.SAU", quote=F, row.names=F, col.names=c("CHROM","POS")) #On ecrit dans le fichier qui sers ensuite à faire le tris dans vcf tools


