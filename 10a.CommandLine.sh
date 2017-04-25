################################################################################################################################
#			Command lines for RAD seq analysis FOR DADI Analysis
################################################################################################################################
#I suggest carefull reading of the lines and modification according to your own needs before running it in a #!/bin/bash

libR=#path to the folder containing the Rscripts

################################################################################################################################
#Renommer et changer les vcf de place#
list_pop=list_pop  #a file containing the list of all pop with one pop by row
for i in $(cat list_pop); do 
    cp $i/batch_1.sumstats.tsv batch_1.sumstats.$i.tsv 
    cp $i/batch_1.vcf batch_1.$i.vcf ; 
done

#Preparation des listes d'individus par pop
for i in $(cat list_pop) ; do 
	awk '{print $1 }'      PopMap.$i.Lib3 >> LF_LP.$i 
	awk '$2==1 {print $1}' PopMap.$i.Lib3 >> LF.$i 
	awk '$2==2 {print $1}' PopMap.$i.Lib3 >> LP.$i ;
done

# Récupérer toutes les lignes du fichier sumstat, moins celles d'en tête commençant par un #, pour récupérer les columns  Couper les 5 premières colonnes 
# Les lignes sont répétées 2-3 fois car il y a 2(3) pops donc on ne récupère que les lignes uniques puis la colonne des positions que l'on va coller au fichier vcf 
for i in $(cat list_pop) ; do
	grep  "^#" batch_1.$i.vcf > entete.$i
	grep -v "^#" batch_1.sumstats.$i.tsv | cut -f 1-5 | uniq | cut -f 5 > columns.$i
done

for i in $(cat list_pop); do
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
for i in $(cat list_pop) ; do
	vcftools --vcf batch_1_trim.$i.vcf --keep LF_LP.$i --out batch_1_trim2.$i --recode ;
done

rm batch_1_trim2*
#############################################################################################################################
#Filtre Hardy-Weinberg et %de genotypage
for i in $(cat list_pop); do 
	vcftools --vcf batch_1_trim.$i.vcf --keep LF.$i --geno 0.80  --hwe 0.05 --out batch_1_LF.$i --recode 
	vcftools --vcf batch_1_trim.$i.vcf --keep LP.$i --geno 0.80  --hwe 0.05 --out batch_1_LP.$i --recode ;
done

mkdir LOG 
mv *log LOG/

mkdir BATCH1
mv  batch_1_LF* batch_1_LP* BATCH1/

cd BATCH1

#Récuperer les positions des SNPs 
libR=/path_To_Rscripts
$libR/10b.GrepPolymorphism.R

########## Filtre sur le vcf d'intêter
mv ../batch_1_trim.*.vcf .
cp ../pop_list .
for i in $(cat list_pop)  ; do 
	vcftools --vcf batch_1_trim.$i.vcf --positions pos_commun.$i --out batch_1.$i.commun --recode 
done

mkdir LOG
mv *log *idx LOG/
##############################################################################################################################
#                         Filtre sur MAF       
for i in $(cat list_pop) ; do 
	vcftools --vcf batch_1_LF.$i.recode.vcf --freq --out lf.$i.freq 
	vcftools --vcf batch_1_LP.$i.recode.vcf --freq --out lp.$i.freq ;
	cut  -f2-6 lf.$i.freq.frq | grep -v POS | sed -re 's/:/\t/g' > $i.lf.freq
	cut  -f2-6 lp.$i.freq.frq | grep -v POS | sed -re 's/:/\t/g' > $i.lp.freq ;  
done

#calcul des fréquences sur le fichier global:
for i in $(cat list_pop) ; do
  	vcftools --vcf  batch_1.$i.commun.recode.vcf --freq --out $i.global.freq
	cut  -f2-6 $i.global.freq.frq > $i.glob.freqglob
	grep -v POS $i.glob.freqglob | sed -re 's/:/\t/g' > $i.glob.freqglob1
done

mv *log *idx LOG/

$libR/10c.MAF.R

mkdir FREQ
mv freql* *.frq FREQ
#CREATION DU DATASET Filtré sur MAF, HWE et %ind génotypés & Calcul de la mean depth
for i in $(cat list_pop) ; do
	vcftools --vcf batch_1_trim.$i.vcf  --positions $i.position_glob.filtered --out batch_LF_LP.$i --recode 
	vcftools --vcf batch_LF_LP.$i.recode.vcf --site-mean-depth --out $i.coverageGLOB_FILTER #mean depth
done

#plot data in R: On regarde la distribution des couvertures
$libR/10d.readDepth.R
############################################################################################################################
#filtre sur la couverture (10à 100)##DATASET DEFINITIF
for i in $(cat list_pop) ; do
	vcftools --vcf batch_LF_LP.$i.recode.vcf --min-meanDP 10 --max-meanDP 100 --out dataLFLP.$i.DEF --recode 
done
#Some verifications in R####################################################################################################
mv ../*sumstat* .
$libR/10e.QualityChecks.R
############################################################################################################################
#Remove heteorzygote excess
$libR/10f.Filtre_HETERO.R #Should use the other script instead (using pegas)

#With pegas:
#first split the dataset
for i in $(cat list_pop) ; do mv ../LF.$i ../LP.$i . ; done
for i in $(cat list_pop) ; do
	vcftools --vcf dataLFLP.$i.DEF.recode.vcf --keep LF.$i --out data.LF.DEF.$i  --recode
	vcftools --vcf dataLFLP.$i.DEF.recode.vcf --keep LP.$i --out data.LP.DEF.$i --recode
done
$libR/10f.filtre_hetero_pegas.R
###########################################################################################################################
#Creates here final dataset
#Using  script Filter_HETERO.R
for i in $(cat list_pop) ; do # CEN JAL ODO OIR RIS SAU ; do
	vcftools --vcf dataLFLP.$i.DEF.recode.vcf  --positions $i.pos_het --out dataLFLP.$i.DEF2 --recode 
done

mv *idx *log LOG

#Using script filter_hetero_pegas.R
#for i in  ; do # CEN JAL ODO OIR RIS SAU ; do
#	vcftools --vcf dataLFLP.$i.DEF.recode.vcf  --positions $i.pos_def.PEGAS --out dataLFLP.$i.DEF2 --recode 
#done

#mv *idx *log LOG
# Creating separates LF & LP final vcf ##################################################################################
for i in $(cat list_pop) ; do cp ../LF.$i ../LP.$i . ; done

for i in $(cat list_pop)  ; do
	vcftools --vcf dataLFLP.$i.DEF2.recode.vcf --keep LF.$i --out data.LF.DEF2.$i  --recode
	vcftools --vcf dataLFLP.$i.DEF2.recode.vcf --keep LP.$i --out data.LP.DEF2.$i --recode
done

mv *idx *log LOG
##########################################################################################################################
#Puis calculer les fréquences alléliques pour chaque VCF,
for i in $(cat list_pop) ; do
	vcftools --vcf data.LF.DEF2.$i.recode.vcf  --freq --out dataset_LF_freq.$i
	vcftools --vcf data.LP.DEF2.$i.recode.vcf --freq --out dataset_LP_freq.$i
done
#Calculer la diversité nucléotidique par site ##
for i in $(cat list_pop) ; do
	vcftools --vcf data.LF.DEF2.$i.recode.vcf --site-pi --out Pi_LF.$i
	vcftools --vcf data.LP.DEF2.$i.recode.vcf --site-pi --out Pi_LP.$i
done
# Computing Fst
for i in $(cat list_pop)  ; do 
    for j in $(cat list_pop) ; do 
        vcftools --vcf dataLFLP.$i.DEF2.recode.vcf --weir-fst-pop $i --weir-fst-pop $j --out fst.$i.vs.$j ; 
    done ; 
done

for i in *log ; do 
    grep "weighted"  $i |sed -e 's/Weir and Cockerham weighted Fst estimate: //g' >> table.$i ; 
done

cat table.fst.* > table.fst
ls fst.pop.*log > list
paste list table.fst > table.fst.pop
sed -i -e 's/fst.pop.//g'  -e 's/log//g' -e 's/pop//g' table.fst.pop -e  's/.\t/\t/g' table.fst.pop 
