#!/bin/bash

#script to Run ustacks
#WARNING: this contains two version, comment out the one that you do not wish to use
# recuperation de l'email
if [ ! -d "03-stacks/" ]  ; then  
    echo "creation du dossier" ; 
    mkdir 03-stacks/; 
fi
mkdir 02-demultiplex/clone_filter/R1
mkdir 02-demultiplex/clone_filter/R2
mv 02-demutliplex/clone_filter/*R1* 02-demultiplex/clone_filter/R1
mv 02-demutliplex/clone_filter/*R2* 02-demultiplex/clone_filter/R2

EMAIL="yourmail@yourmail"
QUEUE=$1
shift 1 

##########################################################################################################################
cpt=1
for FILE in $*; do
    echo "ustacks -t fastq -f $FILE -i $cpt -o 03-stakcs/ -m 4 -M 4 -N 4 -d -r -p 8 --max_locus_stacks 4" >> ustacks.qarray 
    let cpt+=1
done

qarray -m bea -V -cwd -M $EMAIL -e /log/ustacks -q $QUEUE  ustacks.qarray

#see : http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php for more options to passed to ustacks
##########################################################################################################################
#Alternatively using a qsub instead of qarray (if not available on your cluster)
# generation du fichier qsub
mkdir 03-stacks

echo "#!/bin/bash" > ustacks.qsub
echo "#$ -S /bin/bash" >> ustacks.qsub
echo "#$ -m bea" >> ustacks.qsub
echo "#$ -V" >> ustacks.qsub
echo "#$ -cwd" >> ustacks.qsub
echo "#$ -M $EMAIL" >> ustacks.qsub
echo "#$ -e /log/ustacks" >> ustacks.qsub
echo "outfile=03-stacks" >> ustacks.qsub
echo "ustacks -t fastq -f \$INPUT -i \$ID -o $outfile -m 4 -M 3 -N 5 -d -r -p 8 --max_locus_stacks 4" >> ustacks.qsub

# boucle sur qsub
cpt=1
for FILE in 02-demultiplex/clone_filter do
    qsub -o /home/qrougemont/work/RAD/Stacks/log/ustacks -q $QUEUE -v INPUT=$FILE -v ID=$cpt ustacks.qsub 
    let cpt+=1
done
##########################################################################################################################
