#!/bin/bash

for f in lib2 lib3 lib4 lib5 lib6 lib7 lib8 ; do java -jar /home/ubuntu/RAD_Lamproies/Analyses/GBSX-master/releases/latest/GBSX_v1.0.1.jar --Demultiplexer -f1 /home/ubuntu/RAD_Lamproies/Analyses/RawDataMGX/$f.R1.fastq.gz  -f2 /home/ubuntu/RAD_Lamproies/Analyses/RawDataMGX/$f.R2.fastq.gz -i /home/ubuntu/RAD_Lamproies/Analyses/barcodes/barcodes.$f.txt -gzip true  -rad true -o /home/ubuntu/RAD_Lamproies/Analyses/demultiplex/$f; done
done

#Edit: runs faster in parrallel mode than in a loop.
#For instance using 8 qsub on a cluster...

