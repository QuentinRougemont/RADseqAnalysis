################################################################################################################################
####			Command lines for RAD seq analysis FOR DADI Analysis
####
################################################################################################################################
###I suggest carefull reading of the lines and modification according to your own needs before running it in a #!/bin/bash

libR=#path to the folder containing the Rscripts

################################################################################################################################
#Renommer et changer les vcf de place#
for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do cp $i/batch_1.sumstats.tsv batch_1.sumstats.$i.tsv ; done
for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do cp $i/batch_1.vcf batch_1.$i.vcf ; done

#Preparation des listes d'individus par pop
for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do 
	awk '{print $1 }' PopMap.$i.Lib3 >> LF_LP.$i 
	awk '$2==1 {print $1}' PopMap.$i.Lib3 >> LF.$i 
	awk '$2==2 {print $1}' PopMap.$i.Lib3 >> LP.$i ;
done

## Récupérer toutes les lignes du fichier sumstat, moins celles d'en tête commençant par un #, pour récupérer les columns  Couper les 5 premières colonnes 
## Les lignes sont répétées 2-3 fois car il y a 2(3) pops donc on ne récupère que les lignes uniques puis la colonne des positions que l'on va coller au fichier vcf 
for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do
	grep  "^#" batch_1.$i.vcf > entete.$i
	grep -v "^#" batch_1.sumstats.$i.tsv | cut -f 1-5 | uniq | cut -f 5 > columns.$i
done


for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do
	paste <(awk '{print $0}' columns.$i) <(grep -v "^#" batch_1.$i.vcf) > batch_1_tagged.$i.vcf ## Récupérer les lignes du fichier vcf, moins celles d'en tête commençant par un #, et coller columns en première col ##
	awk '$1>5 {print $0}' batch_1_tagged.$i.vcf > batch_1_6.$i.vcf ## Virer pos 1 à 6
	awk '$1!=59 {print $0}' batch_1_6.$i.vcf > batch_1_6_59.$i.vcf ## virer position 59 qui a merdé dans la lib3
	awk '$1<81 {print $0}' batch_1_6_59.$i.vcf > batch_1_6_59_80.$i.vcf ## virer les position > 80
	cut -f 2-150 batch_1_6_59_80.$i.vcf > batch_1_trimmed.$i.vcf ## Effacer la première colonne du vcf pour revenir au format de base ##
	cat  entete.$i batch_1_trimmed.$i.vcf >> batch_1_trim.$i.vcf ## Recoler l'en-tête des vcf , c'est terminé : ##
done
rm entete* *tagged* batch_1_6* batch_1_trimmed*


############################################################################################################################
#Enlever les lamproies marines des jeu de données si elles ont été exportées à l'aide du module stacks:
for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do
	vcftools --vcf batch_1_trim.$i.vcf --keep LF_LP.$i --out batch_1_trim2.$i --recode ;
done
for i in BRE RIS; do mv batch_1_trim2.$i.recode.vcf batch_1_trim.$i.vcf ; done

rm batch_1_trim2*
###############################################################################################################################



#Filtre profondeur #A faire plus loin####
#for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do 
#	vcftools --vcf batch_1_trim.$i.vcf --min-meanDP 10 --max-meanDP 100 --out batch_1_trim0.$i --recode 
#done



#Filtre Hardy-Weinberg et %de genotypage
for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do 
	vcftools --vcf batch_1_trim.$i.vcf --keep LF.$i --geno 0.80  --hwe 0.05 --out batch_1_LF.$i --recode 
	vcftools --vcf batch_1_trim.$i.vcf --keep LP.$i --geno 0.80  --hwe 0.05 --out batch_1_LP.$i --recode ;
done

mkdir LOG 
mv *log LOG/

mkdir BATCH1
mv  batch_1_LF* batch_1_LP* BATCH1/


cd BATCH1


#########Récuperer les positions des SNPs 
libR=/home/ubuntu/Desktop/Script/Rscripts
$libR/10b.GrepPolymorphism.R


########## Filtre sur le vcf d'intêter
mv ../batch_1_trim.*.vcf .

for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do 
	vcftools --vcf batch_1_trim.$i.vcf --positions pos_commun.$i --out batch_1.$i.commun --recode 
done

mkdir LOG
mv *log *idx LOG/

########################################################
##                         Filtre sur MAF       
########################################################
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do 
	vcftools --vcf batch_1_LF.$i.recode.vcf --freq --out lf.$i.freq 
	vcftools --vcf batch_1_LP.$i.recode.vcf --freq --out lp.$i.freq ;
	cp lf.$i.freq.frq  freqlf.$i 
	cp lp.$i.freq.frq  freqlp.$i 
	cut  -f2-6 freqlf.$i | grep -v POS | sed -re 's/:/\t/g' > $i.lf.freq
	cut  -f2-6 freqlp.$i | grep -v POS | sed -re 's/:/\t/g' > $i.lp.freq ;  
done

#calcul des fréquences sur le fichier global:
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
  	vcftools --vcf  batch_1.$i.commun.recode.vcf --freq --out $i.global.freq
	cut  -f2-6 $i.global.freq.frq > $i.glob.freqglob
	grep -v POS $i.glob.freqglob | sed -re 's/:/\t/g' > $i.glob.freqglob1
done

mv *log *idx LOG/

$libR/10c.MAF.R

mkdir FREQ
mv freql* *.frq FREQ


#CREATION DU DATASET Filtré sur MAF, HWE et %ind génotypés & Calcul de la mean depth
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf batch_1_trim.$i.vcf  --positions $i.position_glob.filtered --out batch_LF_LP.$i --recode 
	vcftools --vcf batch_LF_LP.$i.recode.vcf --site-mean-depth --out $i.coverageGLOB_FILTER #mean depth
done


#plot data in R: On regarde la distribution des couvertures

$libR/10d.readDepth.R


######################################################################
#filtre sur la couverture (5 à 60) ##DATASET DEFINITIF
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf batch_LF_LP.$i.recode.vcf --min-meanDP 10 --max-meanDP 100 --out dataLFLP.$i.DEF --recode 
done


##################Some verifications in R##############################
mv ../*sumstat* .
$libR/10e.QualityChecks.R


#########################################################################

#Remove heteorzygote excess
$libR/10f.Filtre_HETERO.R #Should use the other script instead (using pegas)


#With pegas:
#first split the dataset
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do mv ../LF.$i ../LP.$i . ; done
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf dataLFLP.$i.DEF.recode.vcf --keep LF.$i --out data.LF.DEF.$i  --recode
	vcftools --vcf dataLFLP.$i.DEF.recode.vcf --keep LP.$i --out data.LP.DEF.$i --recode
done


$libR/10f.filtre_hetero_pegas.R





############################################################################
#############Creates here final dataset######################################

#Using  script Filter_HETERO.R
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do # CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf dataLFLP.$i.DEF.recode.vcf  --positions $i.pos_het --out dataLFLP.$i.DEF2 --recode 
done

mv *idx *log LOG

#Using script filter_hetero_pegas.R
#for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do # CEN JAL ODO OIR RIS SAU ; do
#	vcftools --vcf dataLFLP.$i.DEF.recode.vcf  --positions $i.pos_def.PEGAS --out dataLFLP.$i.DEF2 --recode 
#done

#mv *idx *log LOG

################ Creating separates LF & LP final vcf ##################################
for i in AA BET BRE CEN JAL ODO OIR RIS SAU; do cp ../LF.$i ../LP.$i . ; done



for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf dataLFLP.$i.DEF2.recode.vcf --keep LF.$i --out data.LF.DEF2.$i  --recode
	vcftools --vcf dataLFLP.$i.DEF2.recode.vcf --keep LP.$i --out data.LP.DEF2.$i --recode
done

mv *idx *log LOG



########### Computing Fst ###############################
for  i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf dataLFLP.$i.DEF2.recode.vcf --weir-fst-pop LF.$i --weir-fst-pop LP.$i --out LFvsLP.$i --recode
done

mv *idx *log LOG

###############################################################################################
#Puis calculer les fréquences alléliques pour chaque VCF,
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf data.LF.DEF2.$i.recode.vcf  --freq --out dataset_LF_freq.$i
	vcftools --vcf data.LP.DEF2.$i.recode.vcf --freq --out dataset_LP_freq.$i
done


#Calculer la diversité nucléotidique par site ##
for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf data.LF.DEF2.$i.recode.vcf --site-pi --out Pi_LF.$i
	vcftools --vcf data.LP.DEF2.$i.recode.vcf --site-pi --out Pi_LP.$i
done


#Weir-fst-globaldataset
#for i in AA BET BRE CEN JAL ODO OIR RIS SAU ; do
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.AA --weir-fst-pop LP.$i --out LF.BETvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.BET --weir-fst-pop LP.$i --out LF.BREvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.BRE --weir-fst-pop LP.$i --out LF.JALvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.CEN --weir-fst-pop LP.$i --out LF.ODOvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.JAL --weir-fst-pop LP.$i --out LF.vsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.ODO --weir-fst-pop LP.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.OIR --weir-fst-pop LP.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.RIS --weir-fst-pop LP.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.SAU --weir-fst-pop LP.$i --out LF.AAvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.AA --weir-fst-pop  LF.$i --out LF.BETvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.BET --weir-fst-pop LF.$i --out LF.BREvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.BRE --weir-fst-pop LF.$i --out LF.JALvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.CEN --weir-fst-pop LF.$i --out LF.ODOvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.JAL --weir-fst-pop LF.$i --out LF.vsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.ODO --weir-fst-pop LF.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.OIR --weir-fst-pop LF.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.RIS --weir-fst-pop LF.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LF.SAU --weir-fst-pop LF.$i --out LF.AAvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.AA --weir-fst-pop LP.$i --out LF.BETvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.BET --weir-fst-pop LP.$i --out LF.BREvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.BRE --weir-fst-pop LP.$i --out LF.JALvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.CEN --weir-fst-pop LP.$i --out LF.ODOvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.JAL --weir-fst-pop LP.$i --out LF.vsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.ODO --weir-fst-pop LP.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.OIR --weir-fst-pop LP.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.RIS --weir-fst-pop LP.$i --out LFvsLP.$i &>> fst_TOT
#	vcftools --vcf batch_1.FULL.TEST7.recode.vcf --weir-fst-pop LP.SAU --weir-fst-pop LP.$i --out LF.AAvsLP.$i &>> fst_TOT
#done

#########################################################################################################################################
#Edit: For population were hybrids were identified, it is important to create new vcf for subsequent analysis!
#They have to be removed for dadi inference.
#If one wants to perform some inference with dadi you MUST NOT used a MAF! the filtering process is different. I will update it one day...

##########################################################################################################################################
#						Example analysis
#						GENET POP
#
###########################################################################################################################################
#TO DO

