#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=salmon_quant
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mili7948@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem=20gb
#SBATCH --time=04:30:00
#SBATCH --output=logs/salmon_fixed_%A_%a.out
#SBATCH --error=logs/salmon_fixed_%A_%a.err
#SBATCH --array=1-638%50

echo "=== Salmon quant ==="
date
echo "ID: $SLURM_ARRAY_TASK_ID"

export SINGULARITY_CACHEDIR=/scratch/Shares/rinn/ML/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR

CSV="/scratch/Shares/rinn/ML/hipsc_timecourse_GFP/rawdata/samples.csv"
OUTDIR="/scratch/Shares/rinn/ML/salmon_with_bootstraps/salmon_quant"
THREADS=7

mkdir -p logs
mkdir -p "$OUTDIR"


LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${LINE_NUM}p" "$CSV")
SAMPLE=$(echo "$LINE" | cut -d',' -f1)
R1=$(echo "$LINE" | cut -d',' -f2)
R2=$(echo "$LINE" | cut -d',' -f3)

echo "sample: $SAMPLE"
echo "Raw R1: $R1"
echo "Raw R2: $R2"

SAMPLE_OUT="$OUTDIR/${SAMPLE}_quant"

if [[ -f "$SAMPLE_OUT/quant.sf" ]]; then
    echo "output exist, skip: $SAMPLE"
    exit 0
fi


R1_CONTAINER="/data$(echo "$R1" | sed 's|/scratch/Shares/rinn/ML/hipsc_timecourse_GFP||')"
R2_CONTAINER="/data$(echo "$R2" | sed 's|/scratch/Shares/rinn/ML/hipsc_timecourse_GFP||')"

echo "vessel R1: $R1_CONTAINER"
echo "vessel R2: $R2_CONTAINER"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "error: file not exist"
    exit 1
fi

mkdir -p "$SAMPLE_OUT"
echo "output directory: $SAMPLE_OUT"
echo "threads: $THREADS"


echo "=== Run Salmon ==="

singularity exec \
    --bind /scratch/Shares/rinn/ML/hipsc_timecourse_GFP:/data \
    --bind /scratch/Shares/rinn/ML/salmon_with_bootstraps:/output \
    docker://combinelab/salmon:latest \
    salmon quant \
    --libType A \
    --index /data/genomes/salmon_index_v43 \
    --mates1 "$R1_CONTAINER" \
    --mates2 "$R2_CONTAINER" \
    --output "/output/salmon_quant/${SAMPLE}_quant" \
    --threads "$THREADS" \
    --dumpEq \
    --posBias --seqBias --gcBias \
    --numGibbsSamples 100

EXIT_CODE=$?
echo "exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ success: $SAMPLE"
    if [ -f "$SAMPLE_OUT/quant.sf" ]; then
        echo "✓ output file created"
        echo "file preview:"
        head -3 "$SAMPLE_OUT/quant.sf"
    fi
else
    echo "✗ fail: $SAMPLE"
    exit 1
fi

echo "=== Done ==="
date
