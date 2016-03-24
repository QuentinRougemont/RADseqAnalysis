#!/bin/bash

for file in lib1 lib2 lib3 lib4 lib5 lib6 lib7 lib8 ; do
	ls $file/*.1.fq.gz > $file/listR1 &&
	ls $file/*.2.fq.gz > $file/listR2 &&
	src=/home/ubuntu/RAD_Lamproies/Analyses/MGX/02a_Demultiplex/ &&
	sro=/home/ubuntu/RAD_Lamproies/Analyses/MGX/02b_CloneFilter/ &&
	while read R1 ;  do echo "clone_filter -1 $src$R1"; done <$file/listR1 > $file/tmp1 &&
	while read R2 ;  do echo "-2 $src$R2 -o $sro  -i gzfastq  &>> /home/ubuntu/RAD_Lamproies/Analyses/MGX/02b_CloneFilter/log_stacks/clone_filt.log"; done <$file/listR2 > $file/tmp2 &&
	paste -d " "  $file/tmp1 $file/tmp2 >> $file/clone.tmp.sh &&
	echo '#!/bin/bash' > $f/header && #be careful to use single quote, otherwise bash search its own history for !/bin/bash
	cat $f/header $f/Clone.tmp.sh > cloneFilter1.sh &&
	chmod 777 $f/cloneFilter1.sh ;
	#./home/ubuntu/RAD_Lamproies/Analyses/MGX/02a_Demultiplex/cloneFilter1.sh &
done



rm *tmp*
