#!/bin/bash
#SBATCH -J atac_plot
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -o %j.out
#SBATCH -e %j.err
# ==========================================
# 0. Setup
# ==========================================
conda init
conda activate deeptools  

# ==========================================
# 1. Path
# ==========================================

BED_FILE="/scratch/Shares/rinn/ML/ATAC_Footprint/TOBIAS_Results/CDRS_209_Ensembl_Fixed.bed"
BW_CONTROL="/scratch/Shares/rinn/ML/ATAC_Footprint/TOBIAS_Results/TOBIAS_0h/merged_0h_corrected.bw"
BW_TREAT="/scratch/Shares/rinn/ML/ATAC_Footprint/TOBIAS_Results/TOBIAS_2.5h/merged_2.5h_corrected.bw"
OUT_PREFIX="CDRS_209"
THREADS=8

# ==========================================
# 2. computeMatrix
# ==========================================

computeMatrix reference-point \
    --referencePoint TSS \
    -b 2000 -a 2000 \
    -R $BED_FILE \
    -S $BW_CONTROL $BW_TREAT \
    --skipZeros \
    --missingDataAsZero \
    -o matrix_${OUT_PREFIX}.gz \
    --numberOfProcessors $THREADS

# ==========================================
# 3. plotProfile
# ==========================================

plotProfile -m matrix_${OUT_PREFIX}.gz \
    -out Fig_TSS_Profile_${OUT_PREFIX}.pdf \
    --perGroup \
    --colors "#377EB8" "#E41A1C" \
    --samplesLabel "0h (Control)" "2.5h (Dox)" \
    --refPointLabel "TSS" \
    --plotTitle "Chromatin Accessibility at CDRS Promoters" \
    --yAxisLabel "Tn5 Cut Signal (Corrected)" \
    --plotHeight 9 --plotWidth 12

# ==========================================
# 4. plotHeatmap
# ==========================================

plotHeatmap -m matrix_${OUT_PREFIX}.gz \
    -out Fig_TSS_Heatmap_${OUT_PREFIX}.pdf \
    --colorMap Blues Reds \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 --zMax 0.1 \
    --kmeans 1


