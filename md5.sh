cd correct

for i in `ls *fastq.gz`; 
do
   gunzip $i
done

for f in `ls *fastq`;
do
   m1=$(md5sum $f)
   cd ..	
   if [ -f $f ];
   then
      m2=$(md5sum $f)
      if [ "$m1" == "$m2" ];
      then
          echo $f 
	  echo "same"
	  mv correct/$f correct/$f.tmp
       else
          echo $f
	  echo $m1
	  echo $m2
       fi
    else
	echo "missing"
	echo $f
    fi 
    cd correct
done



## mdsum of current
## md5 sum of in ..
##compare
