#!/bin/bash

# recuperation de l'email
EMAIL="quentin.rougemont@rennes.inra.fr"
QUEUE=$1
shift 1 #sert à virer l'arguments n°1



##########################################################################################################################
cpt=1
for FILE in $*; do
    echo "ustacks -t fastq -f $FILE -i $cpt -o /home/qrougemont/work/RAD/Stacks/1_ustacks/RES/ -m 4 -M 3 -N 5 -d -r -p 8 --max_locus_stacks 4" >> ustacks.qarray 
    let cpt+=1
done

qarray -m bea -V -cwd -M $EMAIL -e /home/qrougemont/work/RAD/Stacks/log/ustacks -q $QUEUE  ustacks.qarray

#see : http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php for more options to passed to ustacks

##########################################################################################################################
#Alternative using a qsub instead of qarray (if not available on your cluster
# generation du fichier qsub
echo "#!/bin/bash" > ustacks.qsub
echo "#$ -S /bin/bash" >> ustacks.qsub
echo "#$ -m bea" >> ustacks.qsub
echo "#$ -V" >> ustacks.qsub
echo "#$ -cwd" >> ustacks.qsub
echo "#$ -M $EMAIL" >> ustacks.qsub
echo "#$ -e /home/qrougemont/work/RAD/Stacks/log/ustacks" >> ustacks.qsub
echo "ustacks -t fastq -f \$INPUT -i \$ID -o /home/qrougemont/work/RAD/Stacks/1_ustacks/RES -m 4 -M 3 -N 5 -d -r -p 8 --max_locus_stacks 4" >> ustacks.qsub

# boucle sur qsub
cpt=1
for FILE in $*; do
    qsub -o /home/qrougemont/work/RAD/Stacks/log/ustacks -q $QUEUE -v INPUT=$FILE -v ID=$cpt ustacks.qsub 
    let cpt+=1
done

##########################################################################################################################

