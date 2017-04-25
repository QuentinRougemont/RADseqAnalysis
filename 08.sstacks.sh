#!/bin/bash

##########################################################################################################################
####################### A version I used on my own laptop (cluster) machine (8 cores, 32 gb of RAM) #######################
p="-p 8"                  
                         
b="-b 1"                  
c="-c 03-stacks"  
o="-o 03-stacks"          

# Launch sstacks on all samples
for file in $(ls -1 03-stacks/*.tags.tsv | grep -v catalog | perl -pe 's/\.tags\.tsv//')
do
    sstacks $p $b $c  $o -s $file
done &>> sstacks.log

# see http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php for more options to passed to sstacks
##########################################################################################################################
####################### A cluster version I used on genotoul bionformatic platform                #######################

EMAIL="yourmail@yourmail"
QUEUE=$1
shift 1 

echo "#!/bin/bash" > sstacks.qsub
echo "#$ -S /bin/bash" >> sstacks.qsub
echo "#$ -m bea" >> sstacks.qsub
echo "#$ -V" >> sstacks.qsub
echo "#$ -cwd" >> sstacks.qsub
echo "#$ -M $EMAIL" >> sstacks.qsub
echo "#$ -e log/ustacks" >> sstacks.qsub
echo "sstacks -b 1  -o 03-stacks -c  -c 03-stacks -s \$INPUT  -p 8 " >> cstacks.qsub

# boucle sur qsub
cpt=1
for FILE in $*; do
    qsub -o log/cstacks -q $QUEUE -v INPUT=$FILE -v ID=$cpt cstacks.qsub 
    let cpt+=1
done
