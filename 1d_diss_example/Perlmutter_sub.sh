#!/bin/bash
#SBATCH --job-name=dissKBEtest
#SBATCH --account=m5202
#SBATCH -q premium                       # perlmutter regular queue
#SBATCH -C cpu                # perlmutter CPU partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --time=08:00:00                # default (can be overridden via sbatch CLI)
#SBATCH --output=perlmutter_%j.out
#SBATCH --error=perlmutter_%j.err

export SLURM_CPU_BIND="cores"
srun /global/homes/t/tblommel/Libraries/H-NESSi/build/1d_diss_example/Main_diss.x input.inp /global/cfs/projectdirs/m5202/Zeno/IC_25survive_newg/KBE/g=0.1/_N=4_4h.h5
