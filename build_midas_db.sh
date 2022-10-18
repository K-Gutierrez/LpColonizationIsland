. /u/local/Modules/default/init/modules.sh
module load singularity

indir=/u/scratch/r/rwolff/Evolution_Experiment/midas_db/Lactobacillus
mapfile=/u/scratch/r/rwolff/Evolution_Experiment/midas_db/Lactobacillus.mapfile
outdir=/u/scratch/r/rwolff/Evolution_Experiment/midas_db/Lactobacillus

singularity exec $H2_CONTAINER_LOC/MIDAS-1.3.2.sif build_midas_db.py $indir $mapfile $outdir 