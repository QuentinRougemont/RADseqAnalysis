#!/bin/bash

##########################################################################################################################
####################### A version I used on my own laptop (cluster) machine (8 cores, 32 gb of RAM) #######################


p="-p 8"                  
                         
b="-b 1"                  
c="-c /home/ubuntu/RAD_Lamproies/Analyses/GENOTOOL/04a_cstacks_gbsx"  
o="-o /home/ubuntu/RAD_Lamproies/Analyses/04a_cstacks_gbsx"          

# Launch sstacks on all samples
for file in $(ls -1 /home/ubuntu/RAD_Lamproies/Analyses/04a_cstacks_gbsx/*.tags.tsv | grep -v catalog | perl -pe 's/\.tags\.tsv//')
do
    sstacks $p $b $c  $o -s $file
done &>> sstacks.log

# see http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php for more options to passed to sstacks


##########################################################################################################################
####################### A cluster version I used on genotoul bionformatic platform                #######################

EMAIL="quentin.rougemont@rennes.inra.fr"
QUEUE=$1
shift 1 #sert à virer l'arguments n°1

echo "#!/bin/bash" > sstacks.qsub
echo "#$ -S /bin/bash" >> sstacks.qsub
echo "#$ -m bea" >> sstacks.qsub
echo "#$ -V" >> sstacks.qsub
echo "#$ -cwd" >> sstacks.qsub
echo "#$ -M $EMAIL" >> sstacks.qsub
echo "#$ -e /home/qrougemont/work/RAD/Stacks/log/ustacks" >> sstacks.qsub
echo "sstacks -b 1  -o /home/qrougemont/work/RAD/Stacks/2_cstacks/RES -c  -c /home/qrougemont/work/RAD/Stacks/2_cstacks/RES -s \$INPUT  -p 8 " >> cstacks.qsub

# boucle sur qsub
cpt=1
for FILE in $*; do
    qsub -o /home/qrougemont/work/RAD/Stacks/log/cstacks -q $QUEUE -v INPUT=$FILE -v ID=$cpt cstacks.qsub 
    let cpt+=1
done
