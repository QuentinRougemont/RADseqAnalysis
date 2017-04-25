#!/bin/bash

#miniscript to prepare the architecture, move the file and filter them according to some purity criterion
#We work on a total of 8 library, we have eight fastq.gz file for the Read1 and eight fastq.gz for the read2
#TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

#Organise and prepare the file
#prepare directory to store the data
if [ ! -d "01-data/" ]  ; then  
    echo "creation du dossier" ; 
    mkdir 01-data/; 
fi
mv *fastq.gz 01-data/

#Filter reads based on quality filter:
for f in lib* ; do 
	for i in R1 R2 ; do
	zcat 01-data/$f.$i.fastq.gz | grep -A 3 '^@.*[^:]*:N:[^:]*:' |grep -v "^--$" > 01-data/$f.$i.PF.fastq.gz ; 
	done ;
done
