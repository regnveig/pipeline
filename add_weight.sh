#PBS -l walltime=48:00:00
#PBS -q xl230g9q
#PBS -l mem=32gb,ncpus=8


sample=Fir_WT
chr=X
sample_path=/mnt/scratch/ws/aselsukova/202012111301tmp/Rearr_data/done/$sample/aligned


cd /mnt/storage/home/aselsukova/distributives/anaconda_ctale/bin
source activate

cd ~/distributives/scripts/
python add_weight.py $sample_path/$sample.cool $chr $sample_path 

echo got weights

module load jre/1.8.0

cd /mnt/storage/home/aselsukova/distributives/juicer_tools
java -jar juicer_tools.jar addNorm $sample_path/inter.hic $sample_path/vector.txt
# java -jar juicer_tools.jar addNorm $sample_path/inter.hic $sample_path/vector_visualize.txt

echo done
