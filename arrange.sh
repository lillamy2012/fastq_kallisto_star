

for err in `ls *fastq.gz.err`; 
do 
echo $err
  l=$(awk 'END {print NR}' $err)
echo $l 
 if [[ "$l" -eq 0 ]]; then ## no errors
   filename=${err%.*}
   if [ -f "$filename" ]; then  ## gz in folder   
	mv $filename correct/
	rm $err
  elif [ -f correct/"$filename" ];then 
	rm $err
	fi
fi
done
