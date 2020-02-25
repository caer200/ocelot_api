#!/bin/bash
#SBATCH -t 14-00:00:0   				#Time for the job to run 
#SBATCH --job-name=wtuning			   	#Name of the job
#SBATCH -N 1 					#Number of nodes required
#SBATCH -n 24					#Number of cores needed for the job
#SBATCH -p HAS24M128_L
#SBATCH --account=col_cmri235_uksr		#Name of account to run under

# uncomment if necessary
##SBATCH --mail-type ALL				#Send email on start/end
##SBATCH --mail-user smgo234@uky.edu		#Where to send email
##SBATCH --error=SLURM_JOB_%j.err		#Name of error file
##SBATCH --output=SLURM_JOB_%j.out 		#Name of output file

#Module needed for this Gaussian job
module load ccs/gaussian/g16-A.03/g16-haswell   
#module load ccs/gaussian/g16-A.03/g16-sandybridge
echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST " 

# which conda 

#Gaussian Program execution command 
#g16 H8octahedral.com 
# conda init bash
. ~/.bashrc
conda activate ocelot_scratch
# python dummy.py

python trial.py 1b.xyz
