#!/bin/bash

##########################################################################################################################
####################### A version I used on my own laptop (cluster) machine (8 cores, 32 gb of RAM) #######################

infile=03-stacks
outfile=03-stacks
b="-b 1"           # b: MySQL ID of this batch
o="-o "$outfile"/" # o: output path to write results
m="-m 4"           # m: include tags in catalog that match more than one entry
n="-n 3"           # n: number of mismatches allowed between sample tags when
p="-p 8"           # p: enable parallel execution with num_threads threads

s="$(for file in $(ls -1 $infile/*.tags.tsv | perl -pe 's/\.tags\.tsv//'); do echo -s $file; done)"

# Run cstacks
cstacks $b $s $o $g $m $n $p &>> cstacks.log 

# see http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php for more options to passed to cstacks

##########################################################################################################################
####################### A cluster version I used on genotoul bionformatic platform                #######################

EMAIL="yourmal@yourmail"
QUEUE=$1
shift 1 

echo "#!/bin/bash" > cstacks.qsub
echo "#$ -S /bin/bash" >> cstacks.qsub
echo "#$ -m bea" >> cstacks.qsub
echo "#$ -V" >> cstacks.qsub
echo "#$ -cwd" >> cstacks.qsub
echo "#$ -M $EMAIL" >> cstacks.qsub
echo "#$ -e /log/cstacks" >> cstacks.qsub
echo "outfile=03-stacks/" >> cstacks.qsub
echo "cstacks -b 1  -o "$outfile" -m 6 -n 3 -s \$INPUT  -p 8 " >> cstacks.qsub

# boucle sur qsub
cpt=1
for FILE in $*; do
    qsub -o /log/cstacks -q $QUEUE -v INPUT=$FILE -v ID=$cpt cstacks.qsub 
    let cpt+=1
done
