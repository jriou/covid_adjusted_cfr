#!/bin/bash
#SBATCH --job-name='model14B1'
#SBATCH --partition=all
#SBATCH --array=1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=96:00:00

module load Boost/1.66.0-foss-2018a
./model14 sample num_warmup=500 num_samples=500 \
      adapt delta=0.8 algorithm=hmc engine=nuts max_depth=10 init=0.5 \
      data file=data_S_model14B1_2020-03-11-16-34-03.R output file=S_model14B1_2020-03-11-16-34-03_${SLURM_ARRAY_TASK_ID}.csv refresh=10
