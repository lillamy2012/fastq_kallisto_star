
module load sratoolkit/2.8.1-3

for i in `ls *fastq.gz`; 
 do
echo $i 
 gunzip -t $i 2> $i.err
 done
 find . -name "*err" -type f -size +0c -exec -larth {} \;


