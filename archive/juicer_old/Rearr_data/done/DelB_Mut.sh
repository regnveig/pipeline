#PBS -l walltime=48:00:00
#PBS -q xl230g9q
#PBS -l mem=32gb,ncpus=8
TASK=Rearr_data
SAMPLE=DelB_Mut
GENOME=hg38
RESTR=DpnII

HOME=/mnt/storage/home/aselsukova
JU=$HOME/distributives/juicer
TOPDIR=$JU/scratch/$TASK/$SAMPLE
REF=$HOME/genomes/references/$GENOME/$GENOME.fa
RESTRS=$HOME/genomes/restriction_sites/$GENOME"_"$RESTR.txt
CHR=$HOME/genomes/chrom.sizes/chrom.sizes"_"$GENOME
module load jre/1.8.0
cd $JU/scripts
./juicer.sh -t 6 -g $GENOME -D $JU -d $TOPDIR -s $RESTR -p $CHR -y $RESTRS -z $REF 