#PBS -l walltime=52:00:00
#PBS -q bl2x220g7q
#PBS -l mem=24gb,ncpus=12

SAMPLE=DelB
PATH_JSON=path/to/sample/Jsone
OUT_PATH =path/to/output
WDL_PATH+path/to/wdl
echo $SAMPLE

source /mnt/storage/home/aselsukova/miniconda3/bin/activate 
conda activate encode-chip-seq-pipeline
cd ~/scratch/ChIP-Seq_data/$SAMPLE
caper run $WDL_PATH -i $PATH_JSON/$SAMPLE.json --out-dir $OUT_PATH 