. /u/local/Modules/default/init/modules.sh
module load singularity

indir=/u/project/ngarud/rwolff/Evolution_Experiment/midas_db/Lactobacillus
mapfile=/u/project/ngarud/rwolff/Evolution_Experiment/midas_db/Lactobacillus.mapfile
outdir=/u/project/ngarud/rwolff/Evolution_Experiment/midas_db/Lactobacillus

singularity exec $H2_CONTAINER_LOC/MIDAS-1.3.2.sif build_midas_db.py $indir $mapfile $outdir 