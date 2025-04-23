#!/bin/bash

# ----- INPUTS -----
# the directory with preprocessed presence-absence + predictors table
spec_tab_dir='/scratch/pb463/spec_pred_ext/'

# where to save results (will be created if it doesn't exist)
procDir='/scratch/pb463/SDM/'

# Comma separated list of different GEDI scenarios. Otions are 'nogedi', 'wgedicc', 'wgedikr', 'ogedicc', 'ogedikr'
gedi_scens='ogedicc,base,base_lsccdc,base_gedicc,hind_wgedicc'

# Comma separated list of random seeds (corresponding to bootstraps)
seeds='1,2,3,4,5,6,7,8,9,10'

# Comma separated list of feature selection methods
# options 'mrmrcdr,gvifcdr'
fs_methods='gvif,mrmr'

# class balance 
# options: 'bal,nobal'
rf_class_bal='nobal,bal'

# Comma separated list of different regions to model over
# currently not used since site coordinates have been removed
regions='global'
local_ext='117.67,118.15,4.21,4.73'
regional_ext='116.61,118.77,4.15,6.0'
global_ext='80,120,-10,30'



## SLURM settings
# amount of time needed per array job
# 04:00:00 seems ok 
timeReq="02:00:00"

# amount of memory needed per array job
# 4GB seems safe for GLM with dredge
memoryReq="4G"

# number of cores to use
cores=1

# which array jobs to run
# 'all' for all jobs
# comma separated list for select jobs (e.g. '1-10, 15-44' OR '1,3,19,41' OR 'all')
array_jobs_2run='all'
#'all'
#'4401-4500,9601-9700,14801-14900'

# Specify how many array jobs to try and do at once
# 1000
arrayLim=500



# ----- PROCESSING -----
datetime=$(date '+%d%m%Y-%H%M%S')

# Make the output directory structure
if [ ! -d $procDir ]
then
  mkdir -p $procDir
fi

outDir=$procDir's01/'
if [ ! -d $outDir ] 
then
  mkdir -p $outDir/{uni,fs,pred,glm_coefs,eval,imp,pdp,rftrees}
  mkdir -p $procDir's00/'
fi

# Make the slurm output directories
slurmDir=$procDir'slurm/'
if [ ! -d $slurmDir ] 
then
  mkdir -p $slurmDir/{logs,scripts}
fi

# make the run combination file
# column order should be: runid, spec_tabs, seed, method, vif, mrmr, boruta, dredge, parsimony
param_file=$procDir's00/sdm_method_combos.csv'
echo "----- Starting Script 0: building parameter combination table -----"
echo $param_file


runid=1
spec_tabs=$(ls $spec_tab_dir*.csv | tr "\n" ",")
if [ ! -f $param_file ]
then
  touch $param_file
  for r in ${regions//,/ }
  do
    for t in ${spec_tabs//,/ }
    do
        for g in ${gedi_scens//,/ }
        do
            for m in ${fs_methods//,/ }
            do 
                for b in ${rf_class_bal//,/ }
                do
                    for s in ${seeds//,/ }
                    do
                        echo $runid','$r','$t','$g','$m','$b','$s >> $param_file
                        runid=$(( runid + 1 ))
                    done
                done
            done
        done
    done
  done
fi
numFiles=$(wc -l $param_file | cut -d " " -f1)
echo 

# ----- SLURM ARRAY -----

# Create the slurm array job script in the slurm folder
script_s01=$slurmDir'scripts/sbatch_sdm_s01_'$datetime'.sh'

# Specify the basename for the output log file
slurmLog_s01=$slurmDir'logs/sdm_s01_'$datetime'_%A_%a.out'

if [ "$array_jobs_2run" = "all" ]; then
  array_input="1-"$numFiles"%"$arrayLim
else
  array_input=$array_jobs_2run
fi

cat > $script_s01 <<EOT
#!/bin/bash

# - SLURM -
#SBATCH --job-name=sdm_s01
#SBATCH --output=$slurmLog_s01
#SBATCH --time=$timeReq
#SBATCH --mem=$memoryReq
#SBATCH --array=$array_input
# if multiple cores needed in the future set --cpus-per-task=$cores


# get parameters from the param file
# column order should be: runid, region, spec_tab, gedi_scen, fs_meth, seed
p_runid=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $param_file | cut -d "," -f1)
p_reg=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $param_file | cut -d "," -f2)
p_spec_tab=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $param_file | cut -d "," -f3)
p_gedi_scen=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $param_file | cut -d "," -f4)
p_fs_meth=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $param_file | cut -d "," -f5)
p_rfbal=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $param_file | cut -d "," -f6)
p_seed=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $param_file | cut -d "," -f7)

# get the associated region extent
if [ \$p_reg = 'local' ]; then
    p_reg_ext=$local_ext
elif [ \$p_reg = 'regional' ]; then
    p_reg_ext=$regional_ext
elif [ \$p_reg = 'global' ]; then
    p_reg_ext=$global_ext
else
    echo "Invalid name used for region..."
fi

# Activate the environment with the necessary R libraries
module load anaconda3
conda activate /projects/above_gedi/users/pburns/envs/biodiv_mod

# Run the R script
echo "Starting R script..."
echo
srun Rscript /home/pb463/scripts/repos/sunda_sdm_gedi/R/sdm_fs_model_fit_eval.R \$p_runid \$p_reg \$p_reg_ext \$p_spec_tab \$p_gedi_scen \$p_fs_meth \$p_rfbal \$p_seed $cores $outDir 

EOT

# Run the slurm array job script
echo "----- Starting Script 1: sdm_vs_model_fit_eval.R -----"
jobid_s01=$(sbatch --parsable $script_s01)
echo "The SLURM job id is: "$jobid_s01
echo
