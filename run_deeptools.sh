#!/bin/bash
#SBATCH -J atac_plot
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -o %j.out
#SBATCH -e %j.err
# ==========================================
# 0. 环境设置
# ==========================================
conda init
conda activate deeptools  # 修改为你自己的环境名

# ==========================================
# 1. 定义文件路径 (变量化，方便修改)
# ==========================================
# 你的 BED 文件
BED_FILE="/scratch/Shares/rinn/ML/ATAC_Footprint/TOBIAS_Results/CDRS_209_Ensembl_Fixed.bed"

# 你的 bigWig 文件 (请修改为真实文件名)
BW_CONTROL="/scratch/Shares/rinn/ML/ATAC_Footprint/TOBIAS_Results/TOBIAS_0h/merged_0h_corrected.bw"
BW_TREAT="/scratch/Shares/rinn/ML/ATAC_Footprint/TOBIAS_Results/TOBIAS_2.5h/merged_2.5h_corrected.bw"

# 输出文件前缀
OUT_PREFIX="CDRS_209"

# CPU 核心数 (根据服务器情况调整)
THREADS=8

echo "开始运行 deepTools 分析..."
echo "BED: $BED_FILE"
echo "Control: $BW_CONTROL"
echo "Treat: $BW_TREAT"

# ==========================================
# 2. 计算矩阵 (computeMatrix)
# ==========================================
echo "Step 1: Running computeMatrix..."

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
# 3. 画 Profile 图 (plotProfile)
# ==========================================
echo "Step 2: Plotting Profile..."

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
# 4. 画热图 (plotHeatmap)
# ==========================================
echo "Step 3: Plotting Heatmap..."

plotHeatmap -m matrix_${OUT_PREFIX}.gz \
    -out Fig_TSS_Heatmap_${OUT_PREFIX}.pdf \
    --colorMap Blues Reds \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 --zMax 0.1 \
    --kmeans 1

echo "所有分析完成！请下载 PDF 查看结果。"