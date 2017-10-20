module load sratoolkit/2.8.1-3

## in arg = file, column_name, column_type
infile=$1  #"SRARunTable_SCroots.txt"
cn=$2      #15
tn=$3      #7


## single reads  - not here Should be named *_1.fastq

sample_id_single=`awk -v cn="$cn" -v tn="$tn" 'BEGIN {NFR>1} { if ($tn == "SINGLE" || $tn == "SE") print $cn}' $infile`

for sis in $sample_id_single; do
    echo "dumping single read file: " $sis
    fastq-dump $sis
    mv ${sis}.fastq ${sis}_1.fastq
done



## paired
sample_id_paired=`awk -v cn="$cn" -v tn="$tn" 'BEGIN {NFR>1} { if ($tn == "PAIRED" || $tn == "PE") print $cn}' $infile`

for sip in $sample_id_paired; do
   echo "dumping paired read file: " $sip
   fastq-dump --split-files  $sip
done

