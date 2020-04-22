#!/bin/bash
#SBATCH -t 1-00:00:00   				#Time for the job to run 
#SBATCH --job-name=wtuning			   	#Name of the job
#SBATCH -N 1 					#Number of nodes required
#SBATCH -n 16					#Number of cores needed for the job
#SBATCH -p SAN16M64_S
#SBATCH --account=col_cmri235_uksr		#Name of account to run under


#Module needed for this Gaussian job
#module load ccs/gaussian/g16-A.03/g16-haswell   
module load ccs/gaussian/g16-A.03/g16-sandybridge
echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST " 

/scratch/qai222/localpkg/anaconda3/envs/ocelot_fireworks/bin/python trial.py 1b.xyz
