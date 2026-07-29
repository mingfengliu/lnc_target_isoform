#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=hipsc_dox_timecourse
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mili7948@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=36:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "Here we go You've requested $SLURM_CPUS_ON_NODE core."

NXF_SINGULARITY_CACHEDIR=/scratch/Shares/rinn/ZM/singularity_cache
export NXF_SINGULARITY_CACHEDIR

module load singularity/3.1.1

nextflow run nf-core/rnaseq -r 3.14.0 \
-resume \
-profile singularity \
--input /scratch/Shares/rinn/lincxpress/../rnaseq/samplesheet.csv \
--outdir /scratch/Shares/rinn/lincxpress/../rnaseq/nextflow_results/ \
--reads /scratch/Shares/rinn/lincxpress/../rnaseq/fastq/*{_1,_2}.fq.gz \
--fasta /scratch/Shares/rinn/ML/hipsc_timecourse_GFP/genomes/GRCh38.p13.genome.fa \
--gtf /scratch/Shares/rinn/ML/hipsc_timecourse_GFP/genomes/gencode.v38.annotation.gtf \
--pseudo_aligner salmon \
--gencode \
--email mili7948@colorado.edu \
-c nextflow.config

date
