#PBS -q xl230g9q
#PBS -l walltime=1:00:00
#PBS -l mem=32gb,ncpus=8

SAMPLE_NAME=Inv1_Mut
SAMPLE_PATH=~/scratch/Rearr_data/done/$SAMPLE_NAME/aligned/
GENOME=mm10
CHR_SIZES_PATH=/mnt/storage/home/aselsukova/genomes/chrom.sizes/chrom.sizes"_"$GENOME
CHR=chr1
START=70945000
END=81015000


echo $SAMPLE_NAME

source distributives/anaconda_ctale/bin/activate
module load jre/1.8.0

cd $SAMPLE_PATH
awk -F " " '{print $2 "\t" $3 "\t" $6 "\t" $7}' merged_nodups.txt > $SAMPLE_NAME.pairs.txt
cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 $CHR_SIZES_PATH:5000 $SAMPLE_NAME.pairs.txt $SAMPLE_NAME.cool
cp $SAMPLE_NAME.cool $SAMPLE_NAME"_no_balanced".cool
cd /mnt/storage/home/aselsukova/c-tale_normalization/
ctale_normalize $SAMPLE_PATH/$SAMPLE_NAME.cool $CHR":"$START"-"$END output

cd ~/distributives/scripts/
python add_weight.py $SAMPLE_PATH/$SAMPLE_NAME.cool $CHR $SAMPLE_PATH

echo got weights

cd /mnt/storage/home/aselsukova/distributives/juicer_tools
java -jar juicer_tools.jar addNorm $SAMPLE_PATH/inter.hic $SAMPLE_PATH/vector.txt

echo done
