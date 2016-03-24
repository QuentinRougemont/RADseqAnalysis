#!/usr/bin/Rscript

##################################################################################################
#Some lines to perform coplot from Rougemont et al. 2016 Mol Ecol based on results from msstats master


##################################################################################################
#Load data from 
library(extrafont)
#font_import()
loadfonts()


for (i in dir(pattern="_01bis2.eN"))
	{
	tmp=strsplit(i,"_01bis2.eN")[[1]][1]
	tmp2=read.table(i, skip=2)
	assign(tmp, tmp2)
	}	

OIR_en=subset(OIR, OIR[,20]!="NaN")
OIR_en=OIR_en[,20]
max_OIR=max(OIR_en)

AA_en=subset(AA, AA[,20]!="NaN")
AA_en=AA_en[,20]
max_AA=max(AA_en)

BET_en=subset(BET, BET[,20]!="NaN")
BET_en=BET_en[,20]
max_BET=max(BET_en)

RIS_en=subset(RIS, RIS[,20]!="NaN")
RIS_en=RIS_en[,20]
max_RIS=max(RIS_en)

BRE_en=subset(BRE, BRE[,20]!="NaN")
BRE_en=BRE_en[,20]
max_BRE=max(BRE_en)


for (i in dir(pattern="_02bis.eN"))
	{
	tmp=strsplit(i,"_02bis.eN")[[1]][1]
	tmp2=read.table(i, skip=2)
	assign(tmp, tmp2)
	}	

CEN_en=subset(CEN, CEN[,20]!="NaN")
CEN_en=CEN_en[,20]
max_CEN=max(CEN_en)

ODO_en=subset(ODO, ODO[,20]!="NaN")
ODO_en=ODO_en[,20]
max_ODO=max(ODO_en)

SAU_en=subset(SAU, SAU[,20]!="NaN")
SAU_en=SAU_en[,20]
max_SAU=max(SAU_en)

##################################################################################################
#Read Fst data

for(i in dir(pattern="LFvsLP")){
tmp=strsplit(i, ".weir.fst",)[[1]][1] ###à changer probablemenet
tmp2=read.table(i, h=T)
assign(tmp, tmp2)
 }
ls()

LFvsLP.AA= CLEAN.INDEP.LFvsLP.AA
LFvsLP.BET=CLEAN.INDEP.LFvsLP.BET
LFvsLP.OIR=CLEAN.INDEP.LFvsLP.OIR
LFvsLP.RIS=CLEAN.INDEP.LFvsLP.RIS
LFvsLP.BRE=CLEAN.INDEP.LFvsLP.BRE

LFvsLP.AA= CLEAN.LFvsLP.AA
LFvsLP.BET=CLEAN.LFvsLP.BET
LFvsLP.OIR=CLEAN.LFvsLP.OIR
LFvsLP.RIS=CLEAN.LFvsLP.RIS
LFvsLP.BRE=CLEAN.LFvsLP.BRE

LFvsLP.AA[which(LFvsLP.AA[,3]<0),3]<-0
LFvsLP.BET[which(LFvsLP.BET[,3]<0),3]<-0
LFvsLP.OIR[which(LFvsLP.OIR[,3]<0),3]<-0
LFvsLP.RIS[which(LFvsLP.RIS[,3]<0),3]<-0
LFvsLP.RIS[which(LFvsLP.RIS[,3]<0),3]<-0
LFvsLP.BRE[which(LFvsLP.BRE[,3]<0),3]<-0

FstAA_sup90=subset(LFvsLP.AA, LFvsLP.AA[,3] >= max_AA)
FstOIR_sup90=subset(LFvsLP.OIR, LFvsLP.OIR[,3] >= max_OIR)
FstBET_sup90=subset(LFvsLP.BET, LFvsLP.BET[,3] >= max_BET)
FstRIS_sup90=subset(LFvsLP.RIS, LFvsLP.RIS[,3] >= max_RIS)
FstBRE_sup90=subset(LFvsLP.BRE, LFvsLP.BRE[,3] >= max_BRE)


LFvsLP.CEN[which(LFvsLP.CEN[,3]<0),3]<-0
LFvsLP.ODO[which(LFvsLP.ODO[,3]<0),3]<-0
LFvsLP.SAU[which(LFvsLP.SAU[,3]<0),3]<-0

LFvsLP.JAL[which(LFvsLP.JAL[,3]<0),3]<-0

FstCEN_sup90=subset(LFvsLP.CEN, LFvsLP.CEN[,3] >= max_CEN)
FstODO_sup90=subset(LFvsLP.ODO, LFvsLP.ODO[,3] >= max_ODO)
FstSAU_sup90=subset(LFvsLP.SAU, LFvsLP.SAU[,3] >= max_SAU)

#quantJAL.90=quantile(LFvsLP.JAL[,3], c(0.86), na.rm=T) 
#FstJAL_sup90=subset(LFvsLP.JAL, LFvsLP.JAL[,3] >= quantJAL.90) 


A=merge(LFvsLP.OIR,LFvsLP.BET, by="POS")
B=merge(LFvsLP.OIR,LFvsLP.AA , by="POS")
D=merge(LFvsLP.OIR,LFvsLP.RIS, by="POS")
S=merge(LFvsLP.OIR,LFvsLP.BRE, by="POS")
P=merge(LFvsLP.BET,LFvsLP.AA, by="POS")
K=merge(LFvsLP.BET,LFvsLP.RIS, by="POS")
T=merge(LFvsLP.BET,LFvsLP.BRE, by="POS")
R=merge(LFvsLP.AA,LFvsLP.RIS, by="POS")
U=merge(LFvsLP.AA,LFvsLP.BRE, by="POS")
V=merge(LFvsLP.RIS,LFvsLP.BRE, by="POS")


GG=merge(FstOIR_sup90,LFvsLP.BET, by="POS")
KK=merge(FstBET_sup90,LFvsLP.OIR, by="POS")
HH=merge(FstOIR_sup90,LFvsLP.AA, by="POS")
JJ=merge(FstOIR_sup90,LFvsLP.RIS, by="POS")
LL=merge(FstBET_sup90,LFvsLP.AA, by="POS")
MM=merge(FstBET_sup90,LFvsLP.RIS, by="POS")
NN=merge(FstAA_sup90,LFvsLP.OIR, by="POS")
OO=merge(FstAA_sup90,LFvsLP.BET, by="POS")
PP=merge(FstAA_sup90,LFvsLP.RIS, by="POS")
QQ=merge(FstRIS_sup90,LFvsLP.OIR, by="POS")
RR=merge(FstRIS_sup90,LFvsLP.BET, by="POS")
SS=merge(FstRIS_sup90,LFvsLP.AA, by="POS")
TT=merge(FstBRE_sup90,LFvsLP.OIR, by="POS")
TT1=merge(FstOIR_sup90,LFvsLP.BRE, by="POS")
UU=merge(FstBRE_sup90,LFvsLP.AA, by="POS")
UU1=merge(FstAA_sup90,LFvsLP.BRE, by="POS")
VV=merge(FstBRE_sup90,LFvsLP.RIS, by="POS")
VV1=merge(FstRIS_sup90,LFvsLP.BRE, by="POS")
ZZ=merge(FstBRE_sup90,LFvsLP.BET, by="POS")
ZZ1=merge(FstBET_sup90,LFvsLP.BRE, by="POS")

AA1=merge(FstOIR_sup90,FstAA_sup90,by="POS")
BB=merge(FstOIR_sup90,FstBET_sup90,by="POS")
CC=merge(FstOIR_sup90,FstRIS_sup90,by="POS")
DD=merge(FstBET_sup90,FstAA_sup90,by="POS")
EE=merge(FstBET_sup90,FstRIS_sup90,by="POS")
FF=merge(FstAA_sup90,FstRIS_sup90,by="POS")
GGG=merge(FstOIR_sup90,FstBRE_sup90,by="POS")
HHH=merge(FstBET_sup90,FstBRE_sup90,by="POS")
JJJ=merge(FstAA_sup90,FstBRE_sup90,by="POS")
KKK=merge(FstRIS_sup90,FstBRE_sup90,by="POS")


w=lm(BB[,3]~BB[,5])
summary(w)
w1=lm(AA1[,3]~AA1[,5])
summary(w1)
w2=lm(CC[,3]~CC[,5])
summary(w2)
w3=lm(DD[,3]~DD[,5])
summary(w3)
w4=lm(EE[,3]~EE[,5])
summary(w4)
w5=lm(FF[,3]~FF[,5])
summary(w5)
w6=lm(GGG[,3]~GGG[,5])
summary(w6)
w7=lm(HHH[,3]~HHH[,5])
summary(w7)
w8=lm(JJJ[,3]~JJJ[,5])
summary(w8)
w9=lm(KKK[,3]~KKK[,5])
summary(w9)


##################################################################################################
#corrélation

f <- summary(w)$fstatistic
p <- pf(f[1],f[2],f[3],lower.tail=F) #probability function

f1 <- summary(w1)$fstatistic
p1 <- pf(f1[1],f1[2],f1[3],lower.tail=F) #probability function

f2 <- summary(w2)$fstatistic
p2 <- pf(f2[1],f2[2],f2[3],lower.tail=F) #probability function

f3 <- summary(w3)$fstatistic
p3 <- pf(f3[1],f3[2],f3[3],lower.tail=F) #probability function

f4 <- summary(w4)$fstatistic
p4 <- pf(f4[1],f4[2],f3[3],lower.tail=F) #probability function

f5 <- summary(w5)$fstatistic
p5 <- pf(f5[1],f5[2],f5[3],lower.tail=F) #probability function

f6 <- summary(w6)$fstatistic
p6 <- pf(f6[1],f6[2],f6[3],lower.tail=F) #probability function

f7 <- summary(w7)$fstatistic
p7 <- pf(f7[1],f7[2],f7[3],lower.tail=F) #probability function

f8 <- summary(w8)$fstatistic
p8 <- pf(f8[1],f8[2],f8[3],lower.tail=F) #probability function

f9 <- summary(w9)$fstatistic
p9 <- pf(f9[1],f9[2],f9[3],lower.tail=F) #probability function


##################################################################################################
#Parapatry#

K01=merge(LFvsLP.CEN,LFvsLP.JAL, by="POS")
R01=merge(LFvsLP.CEN,LFvsLP.ODO, by="POS")
S01=merge(LFvsLP.CEN,LFvsLP.SAU, by="POS")
T01=merge(LFvsLP.JAL,LFvsLP.ODO, by="POS")
U01=merge(LFvsLP.JAL,LFvsLP.SAU, by="POS")
V01=merge(LFvsLP.ODO,LFvsLP.SAU, by="POS")

W01=merge(LFvsLP.CEN,LFvsLP.BRE, by="POS")
X01=merge(LFvsLP.JAL,LFvsLP.BRE, by="POS")
Y01=merge(LFvsLP.ODO,LFvsLP.BRE, by="POS")
Z01=merge(LFvsLP.SAU,LFvsLP.BRE, by="POS")

EE01=merge(FstCEN_sup90,FstJAL_sup90,by="POS")
FF01=merge(FstCEN_sup90,FstODO_sup90,by="POS")
GG01=merge(FstCEN_sup90,FstSAU_sup90,by="POS")
HH01=merge(FstJAL_sup90,FstODO_sup90,by="POS")
II01=merge(FstJAL_sup90,FstSAU_sup90,by="POS")
KK01=merge(FstODO_sup90,FstSAU_sup90,by="POS")

WW01=merge(FstCEN_sup90,FstBRE_sup90,by="POS")
XX01=merge(FstJAL_sup90,FstBRE_sup90,by="POS")
YY01=merge(FstODO_sup90,FstBRE_sup90,by="POS")
ZZ01=merge(FstSAU_sup90,FstBRE_sup90,by="POS")


CEN2=merge(FstCEN_sup90,LFvsLP.JAL, by="POS")
CEN3=merge(FstCEN_sup90,LFvsLP.ODO, by="POS")
CEN4=merge(FstCEN_sup90,LFvsLP.SAU, by="POS")
JAL2=merge(FstJAL_sup90,LFvsLP.CEN, by="POS")
JAL3=merge(FstJAL_sup90,LFvsLP.ODO, by="POS")
JAL4=merge(FstJAL_sup90,LFvsLP.SAU, by="POS")
ODO2=merge(FstODO_sup90,LFvsLP.CEN, by="POS")
ODO3=merge(FstODO_sup90,LFvsLP.JAL, by="POS")
ODO4=merge(FstODO_sup90,LFvsLP.SAU, by="POS")
SAU2=merge(FstSAU_sup90,LFvsLP.CEN, by="POS")
SAU3=merge(FstSAU_sup90,LFvsLP.JAL, by="POS")
SAU4=merge(FstSAU_sup90,LFvsLP.ODO, by="POS")

BRE20=merge(FstBRE_sup90,LFvsLP.CEN, by="POS")
BRE30=merge(FstBRE_sup90,LFvsLP.JAL, by="POS")
BRE40=merge(FstBRE_sup90,LFvsLP.ODO, by="POS")
BRE50=merge(FstBRE_sup90,LFvsLP.SAU, by="POS")
CEN50=merge(FstCEN_sup90,LFvsLP.BRE, by="POS")
JAL50=merge(FstJAL_sup90,LFvsLP.BRE, by="POS")
ODO50=merge(FstODO_sup90,LFvsLP.BRE, by="POS")
SAU50=merge(FstSAU_sup90,LFvsLP.BRE, by="POS")




w401=lm(EE01[,3]~EE01[,5])
w501=lm(FF01[,3]~FF01[,5])
w601=lm(GG01[,3]~GG01[,5])
w701=lm(HH01[,3]~HH01[,5])
w801=lm(II01[,3]~II01[,5])
w901=lm(KK01[,3]~KK01[,5])

w1001=lm(WW01[,3]~WW01[,5])
w1101=lm(XX01[,3]~XX01[,5])
w1201=lm(YY01[,3]~YY01[,5])
w1301=lm(ZZ01[,3]~ZZ01[,5])

f401 <- summary(w401)$fstatistic
f501 <- summary(w501)$fstatistic
f601 <- summary(w601)$fstatistic
f701 <- summary(w701)$fstatistic
f801 <- summary(w801)$fstatistic
f901 <- summary(w901)$fstatistic

f1001<- summary(w1001)$fstatistic
f1101<- summary(w1101)$fstatistic
f1201<- summary(w1201)$fstatistic
f1301<- summary(w1301)$fstatistic

p401 <- pf(f401[1],f401[2],f401[3],lower.tail=F) #probability function
p501 <- pf(f501[1],f501[2],f501[3],lower.tail=F) #probability function
p601 <- pf(f601[1],f601[2],f601[3],lower.tail=F) #probability function
p701 <- pf(f701[1],f701[2],f701[3],lower.tail=F) #probability function
p801 <- pf(f801[1],f801[2],f801[3],lower.tail=F) #probability function
p901 <- pf(f901[1],f901[2],f901[3],lower.tail=F) #probability function

p1001<-pf(f1001[1],f1001[2],f1001[3],lower.tail=F)
p1101<-pf(f1101[1],f1101[2],f1101[3],lower.tail=F)
p1201<-pf(f1201[1],f1201[2],f1201[3],lower.tail=F)
p1301<-pf(f1301[1],f1301[2],f1301[3],lower.tail=F)
  

#######################################################################################################"
pdf(file="coplot_Parallelism.pdf", 15,15, family="Times New Roman")
#par(bg="grey92")
par(mfrow=c(4,3))
par(mar=c(4,4,2,2),oma=c(1,2,2,0)) 
plot(A[,3]~A[,5], xlab=expression(paste("BET",(italic(F[ST])))), ylab=expression(paste("OIR",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1)) #, pch=19,cex=1)
points(GG[,5],GG[,3], col="burlywood4", pch=19, cex=1)
points(KK[,3],KK[,5], col="black", pch=19, cex=1)
points(BB[,3]~BB[,5], col="blue", pch=19, cex=1) 
abline(w, lty=1, lwd=1, col="blue")
text( x=0, y=0.9, cex=1.3, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w)$r.squared),2)) ) ) )

plot(B[,3]~B[,5], xlab=expression(paste("AA",(italic(F[ST])))), ylab=expression(paste("OIR",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1)) #, pch=19,cex=1)
points(HH[,5],HH[,3], col="burlywood4", pch=19, cex=1)
points(NN[,3],NN[,5], col="black", pch=19, cex=1)
points(AA1[,3]~AA1[,5], col="blue", pch=19, cex=1)
abline(w1, lty=1, lwd=1, col="blue")

text( x=0, y=0.9, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w1)$r.squared),2)) ) ) )

plot(D[,3]~D[,5], xlab=expression(paste("RIS",(italic(F[ST])))), ylab=expression(paste("OIR",(italic(F[ST])))), col="grey", xlim=c(-0.01,1),ylim=c(0,1)) #, pch=19,cex=1)
points(JJ[,5],JJ[,3], col="burlywood4", pch=19, cex=1)
points(QQ[,3],QQ[,5], col="black", pch=19, cex=1)
points(CC[,3]~CC[,5], col="blue", pch=19, cex=1)
abline(w2, lty=1, lwd=1, col="blue")
text( x=0, y=0.9, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w2)$r.squared),2)) ) ) )

plot(P[,3]~P[,5], xlab=expression(paste("AA",(italic(F[ST])))), ylab=expression(paste("BET",(italic(F[ST])))), col="grey", xlim=c(-0.01,1) , ylim=c(0,1)) #, pch=19,cex=1)
points(LL[,5],LL[,3], col="burlywood4", pch=19, cex=1)
points(OO[,3],OO[,5], col="black", pch=19, cex=1)
abline(w3, lty=1, lwd=1, col="blue")
points(DD[,3]~DD[,5], col="blue", pch=19, cex=1)
text( x=0, y=0.9, cex=1.2,pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w3)$r.squared),2)) ) ) )

plot(K[,3]~K[,5], xlab=expression(paste("RIS",(italic(F[ST])))), ylab=expression(paste("BET",(italic(F[ST])))), col="grey", xlim=c(-0.01,1) , ylim=c(0,1)) #, pch=19,cex=1)
points(MM[,5],MM[,3], col="burlywood4", pch=19, cex=1)
points(RR[,3],RR[,5], col="black", pch=19, cex=1)
points(EE[,3]~EE[,5], col="blue", pch=19, cex=1)
abline(w4, lty=1, lwd=1, col="blue")
text( x=0, y=0.9, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w4)$r.squared),2)) ) ) )

plot(R[,3]~R[,5], xlab=expression(paste("RIS",(italic(F[ST])))), ylab=expression(paste("AA",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1)) #, pch=19,cex=1)
points(PP[,5],PP[,3], col="burlywood4", pch=19, cex=1)
points(SS[,3],SS[,5], col="black", pch=19, cex=1)
points(FF[,3]~FF[,5], col="blue", pch=19, cex=1)
abline(w5, lty=1, lwd=1, col="blue")
text( x=0, y=0.9, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w5)$r.squared),2)) ) ) )


plot(R01[,3]~R01[,5], xlab=expression(paste("ODO",(italic(F[ST])))), ylab=expression(paste("LOI",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1))
points(CEN3[,5],CEN3[,3], col="burlywood4", pch=19, cex=1)
points(ODO2[,3],ODO2[,5], col="black", pch=19, cex=1)
points(FF01[,3]~FF01[,5], col="darkgreen", pch=19, cex=1)
abline(w501, lty=1, lwd=1, col="darkgreen")
text( x=0, y=0.95, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w501)$r.squared),2)) ) ) )

plot(S01[,3]~S01[,5], xlab=expression(paste("GAR",(italic(F[ST])))), ylab=expression(paste("LOI",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1))
points(CEN4[,5],CEN4[,3], col="burlywood4", pch=19, cex=1)
points(SAU2[,3],SAU2[,5], col="black", pch=19, cex=1)
points(GG01[,3]~GG01[,5], col="darkgreen", pch=19 ,cex=1)
abline(w601, lty=1, lwd=1, col="darkgreen")
text( x=0, y=0.95, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w601)$r.squared),2)) ) ) )


plot(V01[,3]~V01[,5], xlab=expression(paste("GAR",(italic(F[ST])))), ylab=expression(paste("ODO",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1))
points(ODO4[,5],ODO4[,3], col="burlywood4", pch=19, cex=1)
points(SAU4[,3],SAU4[,5], col="black", pch=19, cex=1)
points(KK01[,3]~KK01[,5], col="darkgreen", pch=19, cex=1)
abline(w901, lty=1, lwd=1, col="darkgreen")
text( x=0, y=0.95, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w901)$r.squared),2) )) ))


plot(W01[,3]~W01[,5], xlab=expression(paste("BRE",(italic(F[ST])))), ylab=expression(paste("LOI",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1)) #, pch=19,cex=1)
points(CEN50[,5],CEN50[,3], col="burlywood4", pch=19, cex=1)
points(BRE20[,3],BRE20[,5], col="black", pch=19, cex=1)
points(WW01[,3]~WW01[,5], col="darkgreen", pch=19, cex=1)
abline(w1001, lty=1, lwd=1, col="darkgreen")
text( x=0, y=0.9, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w1001)$r.squared),2)) ) ) )


plot(Y01[,3]~Y01[,5], xlab=expression(paste("BRE",(italic(F[ST])))), ylab=expression(paste("ODO",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1)) #, pch=19,cex=1)
points(ODO50[,5],ODO50[,3], col="burlywood4", pch=19, cex=1)
points(BRE40[,3],BRE40[,5], col="black", pch=19, cex=1)
points(YY01[,3]~YY01[,5], col="darkgreen", pch=19, cex=1)
abline(w1201, lty=1, lwd=1, col="darkgreen")
text( x=0, y=0.9, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w1201)$r.squared),2)) ) ) )


plot(Z01[,3]~Z01[,5], xlab=expression(paste("BRE",(italic(F[ST])))), ylab=expression(paste("GAR",(italic(F[ST])))), col="grey", xlim=c(-0.01,1), ylim=c(0,1)) #, pch=19,cex=1)
points(SAU50[,5],SAU50[,3], col="burlywood4", pch=19, cex=1)
points(BRE50[,3],BRE50[,5], col="black", pch=19, cex=1)
points(ZZ01[,3]~ZZ01[,5], col="darkgreen", pch=19, cex=1)
abline(w1301, lty=1, lwd=1, col="darkgreen")
text( x=0, y=0.9, cex=1.2, pos=4, as.expression( substitute( r ==r1,
list(r1=round(sqrt(summary(w1301)$r.squared),2)) ) ) )

#mtext("figure 7", side=3)
dev.off()

