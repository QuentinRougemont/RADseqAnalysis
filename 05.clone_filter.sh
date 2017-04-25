#!/bin/bash

#run clone_filter from stacks to remove PCR duplicate (this is important as it may represent a large part of the data)
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

mkdir 02-demultiplex/clone_filter

outfile="02-demultiplex/clone_filter"
infile="02-demultiplex/cutadapt"

ls $infile/*.R1.fq.gz > $infile/listR1 
ls $infile/*.R2.fq.gz > $infile/listR2
src=$infile  
sro=$outfile  

while read R1 ;  do 
    echo "clone_filter -1 "$src"$R1 ";
done <$infile/listR1 > $infile/tmp1 
while read R2 ;  do 
    echo " -2 $src$R2 -o $sro  -i gzfastq  &>> log/log_"$TIMESTAMP"_clone_filter.log "; 
done <$infile/listR2 > $infile/tmp2 

paste -d " "  $infile/tmp1 $infile/tmp2 >> $infile/clone.tmp.sh 
echo '#!/bin/bash' > $infile/header 
cat $infile/header $infile/clone.tmp.sh > clone_filter.sh
chmod +x $infile/cloneFilter1.sh ;

rm $infile/*tmp*
rm $infile/list*
./clone_filter.sh 
