#!/bin/bash

##########################################################################################################################
####################### A version I used on my own laptop (cluster) machine (8 cores, 32 gb of RAM) #######################


b="-b 1"           
P="-P ./04a_cstacks_gbsx"    
M="-M ./PopMap" 
t="-t 5"
r="-r 0.7"
p="-p 2" 
vcf="--vcf"
m="-m 5" 

for f in PopMapAA0 ; do
	populations $b $P $M/$f $m $t $p $m $vcf &>> stacks_populations.log ;
done 

###Repeat this for each pairs of population to get each pairwise populations vcf that will be used for subsequent analysis

# see http://catchenlab.life.illinois.edu/stacks/comp/populations.php for more options to passed to sstacks


