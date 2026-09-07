#!/bin/bash
#SBATCH --job-name=TOBIAS_Footprint
#SBATCH --output=tobias_%j.log
#SBATCH --error=tobias_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48       
#SBATCH --mem=128G               
#SBATCH --time=72:00:00        
#SBATCH --partition=long        

# ==============================================================================
# 1. Setup
# ==============================================================================


source ~/.bashrc  
conda init bash  
source ~/.bashrc 
conda activate tobias_env 


BAM_DIR="/scratch/Shares/rinn/lincxpress/all_atac/results_broad/bwa/merged_library"
OUT_DIR="/scratch/Shares/rinn/ML/ATAC_Footprint/TOBIAS_Results"
GENOME="/scratch/Shares/rinn/ML/hipsc_timecourse_GFP/genomes/GRCh38.p13.genome.fa"
PEAKS="/scratch/Shares/rinn/ML/ATAC_Footprint/consensus_peaks.mLb.clN.bed"
MOTIFS="/scratch/Shares/rinn/ML/ATAC_Footprint/JASPAR2022_CORE_vertebrates_non-redundant.meme"

THREADS=48
mkdir -p ${OUT_DIR}

echo "Starting TOBIAS pipeline at $(date)"

# ==============================================================================
# 1.5 BAM Integrity check
# ==============================================================================
check_bam_integrity() {
    local bam_file=$1
    if samtools quickcheck $bam_file 2>/dev/null; then
        return 0
    else
        return 1
    fi
}

create_valid_bam_list() {
    local pattern=$1
    local output_file=$2
    
    echo "Creating valid BAM list for pattern: $pattern"
    rm -f $output_file
    local valid_count=0
    
    for bam in $(find $BAM_DIR -name "$pattern"); do
        if check_bam_integrity $bam; then
            echo $bam >> $output_file
            ((valid_count++))
        else
            echo "Skipping corrupted file: $(basename $bam)"
        fi
    done
    
    echo "Found $valid_count valid BAM files"
    return $valid_count
}

# ==============================================================================
# 2. BAM Merge
# ==============================================================================
echo ">>> Step 1: Merging VALID BAM files..."

# --- 0h ---
if [ ! -f "${OUT_DIR}/merged_0h.bam" ]; then
    create_valid_bam_list "*_0h_*.bam" "${OUT_DIR}/valid_bam_list_0h.txt"
    count_0h=$(wc -l < ${OUT_DIR}/valid_bam_list_0h.txt 2>/dev/null || echo 0)
    
    if [ $count_0h -gt 0 ]; then
        echo "Merging $count_0h valid BAM files for 0h..."
        samtools merge -@ ${THREADS} -b ${OUT_DIR}/valid_bam_list_0h.txt ${OUT_DIR}/merged_0h.bam
        samtools index -@ ${THREADS} ${OUT_DIR}/merged_0h.bam
        echo "0h BAM merge completed. File size: $(du -h ${OUT_DIR}/merged_0h.bam | cut -f1)"
    else
        echo "ERROR: No valid BAM files found for 0h!"
        exit 1
    fi
else
    echo "Merged 0h BAM already exists. Skipping."
fi

# ---2.5h ---
if [ ! -f "${OUT_DIR}/merged_2.5h.bam" ]; then
    create_valid_bam_list "*_2.5h_*.bam" "${OUT_DIR}/valid_bam_list_2.5h.txt"
    count_25h=$(wc -l < ${OUT_DIR}/valid_bam_list_2.5h.txt 2>/dev/null || echo 0)
    
    if [ $count_25h -gt 0 ]; then
        echo "Merging $count_25h valid BAM files for 2.5h..."
        samtools merge -@ ${THREADS} -b ${OUT_DIR}/valid_bam_list_2.5h.txt ${OUT_DIR}/merged_2.5h.bam
        samtools index -@ ${THREADS} ${OUT_DIR}/merged_2.5h.bam
        echo "2.5h BAM merge completed. File size: $(du -h ${OUT_DIR}/merged_2.5h.bam | cut -f1)"
    else
        echo "ERROR: No valid BAM files found for 2.5h!"
        exit 1
    fi
else
    echo "Merged 2.5h BAM already exists. Skipping."
fi

# ==============================================================================
# 3. Tn5 ATACorrect
# ==============================================================================
echo ">>> Step 2: Running TOBIAS ATACorrect..."

# 0h 
if [ ! -f "${OUT_DIR}/TOBIAS_0h/merged_0h_corrected.bw" ]; then
    echo "Running ATACorrect for 0h..."
    TOBIAS ATACorrect \
        --bam ${OUT_DIR}/merged_0h.bam \
        --genome ${GENOME} \
        --peaks ${PEAKS} \
        --outdir ${OUT_DIR}/TOBIAS_0h \
        --cores ${THREADS} \

else
    echo "0h Correction already done."
fi

# 2.5h
if [ ! -f "${OUT_DIR}/TOBIAS_2.5h/merged_2.5h_corrected.bw" ]; then
    echo "Running ATACorrect for 2.5h..."
    TOBIAS ATACorrect \
        --bam ${OUT_DIR}/merged_2.5h.bam \
        --genome ${GENOME} \
        --peaks ${PEAKS} \
        --outdir ${OUT_DIR}/TOBIAS_2.5h \
        --cores ${THREADS} \
        --split 10 \
        --track_coverage 5000
else
    echo "2.5h Correction already done."
fi

# ==============================================================================
# 4. ScoreBigwig
# ==============================================================================
echo ">>> Step 3: Running TOBIAS ScoreBigwig..."

# 0h 
if [ ! -f "${OUT_DIR}/TOBIAS_0h/merged_0h_footprints.bw" ]; then
    echo "Running ScoreBigwig for 0h..."
    TOBIAS ScoreBigwig \
        --signal ${OUT_DIR}/TOBIAS_0h/merged_0h_corrected.bw \
        --regions ${PEAKS} \
        --output ${OUT_DIR}/TOBIAS_0h/merged_0h_footprints.bw \
        --cores ${THREADS}
else
    echo "0h ScoreBigwig already done."
fi

# 2.5h 
if [ ! -f "${OUT_DIR}/TOBIAS_2.5h/merged_2.5h_footprints.bw" ]; then
    echo "Running ScoreBigwig for 2.5h..."
    TOBIAS ScoreBigwig \
        --signal ${OUT_DIR}/TOBIAS_2.5h/merged_2.5h_corrected.bw \
        --regions ${PEAKS} \
        --output ${OUT_DIR}/TOBIAS_2.5h/merged_2.5h_footprints.bw \
        --cores ${THREADS}
else
    echo "2.5h ScoreBigwig already done."
fi

# ==============================================================================
# 5. BINDetect
# ==============================================================================
echo ">>> Step 4: Running TOBIAS BINDetect (Differential Footprinting)..."

TOBIAS BINDetect \
    --motifs ${MOTIFS} \
    --signals ${OUT_DIR}/TOBIAS_0h/merged_0h_footprints.bw ${OUT_DIR}/TOBIAS_2.5h/merged_2.5h_footprints.bw \
    --genome ${GENOME} \
    --peaks ${PEAKS} \
    --outdir ${OUT_DIR}/BINDetect_Output \
    --cond_names 0h 2.5h \
    --cores ${THREADS} \


echo "Pipeline finished successfully at $(date)"
