#Organise and prepare the file
mkdir lib1 lib2 lib3 lib4 lib5 lib6 lib7 lib8
for f in lib1 lib2 lib3 lib4 lib5 lib6 lib7 lib8; do mkdir $f/Ri $f/R2 ; done
#mv the file in their lib
#EVENTUALLY: creates a single fastq for each lib :
for f in lib1 lib2 lib3 lib4 lib5 lib6 lib7 lib8; do
	for i in R1 R2;  do cat $f/$i/*.gz >> $f.$i.fastq.gz ;
	done ;
done


0. Filter reads based on quality filter:
for f in lib1 lib2 lib3 lib4 lib5 lib6 lib7 lib8 ; do 
	for i in R1 R2 ; do
	zcat $f/$f.$i.fastq.gz | grep -A 3 '^@.*[^:]*:N:[^:]*:' |grep -v "^--$" > $f.$i.PF.fastq.gz ; 
	done ;
done
