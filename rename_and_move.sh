all_fastq=$(find . -name \*.fastq)
for file in $all_fastq
do
fname=$(realpath $file)
echo $fname 
dname=$(realpath $(dirname $file))
mkdir $dname/fastq
mv $fname $dname/fastq/$(awk -F"[_]" '{print $1 "_R" $2 }' <<< $(basename $file))
done