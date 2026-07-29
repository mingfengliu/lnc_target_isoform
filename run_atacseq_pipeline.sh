#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=all_atacseq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mili7948@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=2-00:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core."
echo "Running Nextflow from directory: $(pwd)"

module load openjdk/21.0.1
module load singularity/3.1.1

nextflow run nf-core/atacseq -r 2.1.2 -resume -c nextflow.config --input /scratch/Shares/rinn/lincxpress/all_atac/all_atacseq_design.csv --outdir /scratch/Shares/rinn/lincxpress/all_atac/results_broad --fasta /scratch/Shares/rinn/ML/hipsc_timecourse_GFP/genomes/GRCh38.p13.genome.fa --gtf /scratch/Shares/rinn/ML/hipsc_timecourse_GFP/genomes/gencode.v38.annotation.gtf --blacklist /scratch/Shares/rinn/ML/hipsc_timecourse_GFP/genomes/hg38-blacklist.v2.bed --read_length 150 -profile singularity \   
