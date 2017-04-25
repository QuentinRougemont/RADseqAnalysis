#!/bin/bash

#script to run cutadapt on paired-end rad sequencing data
#remove adatapter and trim sequenced to a desired length
#BEFORE doing that, do not forget to have a look at the quality of the data (use fastqc for instance)
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
NCPU=8 #number of CPU, can be set to any number

mkdir 01-data/cutadapt 01-data/cutadapt_2
mkdir log

input1=$( ls 01-data/*R1.*gz | sed -e 's/01-data\///g' )
#input2=$( ls 01-data/*R2.*gz | sed -e 's/01-data\///g' )
input2="ls -1 01-data/*R2.gz"
output="-o /01-data/cutadapt"
output2="-p /01-data/cutadapt"
adapter="-a AGATCGGAAGAGCG" #name of the first adapter
adapter2="-A AGACCGATCAGAAC" #name of the second adapter 
error="-e 0.1" #error rate 
m="-m 85" #length at which we want to trim the data

if [[ -z "$NCPU" ]]
then
   NCPU=1
fi

ls -1 01-data/*R1
parallel -j $NCPU cutadapt $adapter $adapter2 $output/$input1 $output2/"${input1%.R1.fq.gz}".R2.fq.gz 
        {} $input2 
        $error 
        #$discard_trimmed $untrimmed_output $discard_short
        $m  &>> /log/log_"$TIMESTAMP"_cutadapt_"${i%.fq.gz}".log

#All information on cutadapt is available here: http://cutadapt.readthedocs.org/en/stable/guide.html
#can be cloned from github at: git clone https://github.com/marcelm/cutadapt.git
#for dependencies see: http://cutadapt.readthedocs.org/en/stable/installation.html
