#!/bin/bash

##########################################################################################################################
####################### A version I used on my own laptop (cluster) machine (8 cores, 32 gb of RAM) #######################

b="-b 1"            # b: MySQL ID of this batch
o="-o /home/ubuntu/RAD_Lamproies/Analyses/04a_cstacks_gbsx"    # o: output path to write results

m="-m 6"             # m: include tags in catalog that match more than one entry
n="-n 3"            # n: number of mismatches allowed between sample tags when
                    #   generating the catalog (default 0)
p="-p 8"           # p: enable parallel execution with num_threads threads

s="$(for file in $(ls -1 /home/ubuntu/RAD_Lamproies/Analyses/03a_ustacks_gbsx20_04/*.tags.tsv | perl -pe 's/\.tags\.tsv//'); do echo -s $file; done)"

# Run cstacks
cstacks $b $s $o $g $m $n $p &>> cstacks.log 

# see http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php for more options to passed to cstacks


##########################################################################################################################
####################### A cluster version I used on genotoul bionformatic platform                #######################

EMAIL="quentin.rougemont@rennes.inra.fr"
QUEUE=$1
shift 1 #sert à virer l'arguments n°1

echo "#!/bin/bash" > cstacks.qsub
echo "#$ -S /bin/bash" >> cstacks.qsub
echo "#$ -m bea" >> cstacks.qsub
echo "#$ -V" >> cstacks.qsub
echo "#$ -cwd" >> cstacks.qsub
echo "#$ -M $EMAIL" >> cstacks.qsub
echo "#$ -e /home/qrougemont/work/RAD/Stacks/log/cstacks" >> cstacks.qsub
echo "cstacks -b 1  -o /home/qrougemont/work/RAD/Stacks/2_cstacks/RES -m 6 -n 3 -s \$INPUT  -p 8 " >> cstacks.qsub

# boucle sur qsub
cpt=1
for FILE in $*; do
    qsub -o /home/qrougemont/work/RAD/Stacks/log/cstacks -q $QUEUE -v INPUT=$FILE -v ID=$cpt cstacks.qsub 
    let cpt+=1
done
