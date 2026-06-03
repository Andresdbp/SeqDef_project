#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=seqdef_runtime
#SBATCH --partition=bigmem
#SBATCH --time=08:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=800G
#SBATCH --output=/scratch/user/andresdbp/seqdef/results/runtime_%j.out
#SBATCH --error=/scratch/user/andresdbp/seqdef/results/runtime_%j.err

cd /scratch/user/andresdbp/seqdef
mkdir -p results

module purge
module load gfbf/2024a R/4.4.2
export R_LIBS_USER=/scratch/user/andresdbp/R_libs_seqdef

# single-threaded BLAS for a clean O(n^2) scaling curve (one calculation = one thread)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "Host: $(hostname)  Date: $(date)"
echo "SLURM_MEM_PER_NODE=$SLURM_MEM_PER_NODE  CPUS=$SLURM_CPUS_ON_NODE"
Rscript grace_runtime_compute.R
echo "Exit: $?  Done: $(date)"
