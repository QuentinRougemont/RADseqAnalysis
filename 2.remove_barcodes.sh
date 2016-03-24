#!/bin/bash

input="/home/ubuntu/RAD_Lamproies/Analyses/RawDataMGX/lib1/Librairie1.R1.PF.fastq"
input2="/home/ubuntu/RAD_Lamproies/Analyses/RawDataMGX/lib2/Librairie2.R2.PF.fastq.gz"
output="-o /home/ubuntu/RAD_Lamproies/Analyses/RawDataMGX/lib1/lib1_cut.fastq"
output2="-p /home/ubuntu/RAD_Lamproies/Analyses/RawDataMGX/lib2/lib2.cut.R2.fastq.gz"
adapter="-a AGATCGGAAGAGCG"
Adapter2="-A AGACCGATCAGAAC"
error="-e 0.1"
m="-m 85"

cutadapt $adapter $Adapter2 $output $output2 $input $input2 $error $discard_trimmed $untrimmed_output $m $discard_short &>>  /home/ubuntu/RAD_Lamproies/Analyses/RawDataMGX/log/log_lib2.log


#repeat this for each Library
#All information on cutadapt is available here: http://cutadapt.readthedocs.org/en/stable/guide.html
#can be download in github using: git clone https://github.com/marcelm/cutadapt.git
#for dependencies see: http://cutadapt.readthedocs.org/en/stable/installation.html

#Do not forget to have a look at the quality of the data (use fastqc for instance)

