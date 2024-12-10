{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww16160\viewh11920\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs28 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 #!/bin/bash\
\
# Comandi epigenomics\
\
# Step 1: Download the .bam files (unfiltered) for the 2 replicates and the input from ENCODE and filter them\
# Filter reads with quality score >= 1\
samtools view -bq 1 replicate1_unfiltered.bam > replicate1_filtered_unique.bam\
samtools view -bq 1 replicate2_unfiltered.bam > replicate2_filtered_unique.bam\
samtools view -bq 1 input_unfiltered.bam > input_filtered_unique.bam\
\
# Step 2: Perform quality control using samtools flagstat\
samtools flagstat replicate1_unfiltered.bam\
samtools flagstat replicate1_filtered_unique.bam\
\
samtools flagstat replicate2_unfiltered.bam\
samtools flagstat replicate2_filtered_unique.bam\
\
samtools flagstat input_unfiltered.bam\
samtools flagstat input_filtered_unique.bam\
\
# Step 3: Peak calling with MACS2\
# First, run peak calling for replicate 1 (filtered)\
macs2 callpeak -t replicate1_filtered_unique.bam -c input_filtered_unique.bam -g hs -n replicate1_peakcalling\
\
# Then, run peak calling for replicate 2 (filtered)\
macs2 callpeak -t replicate2_filtered_unique.bam -c input_filtered_unique.bam -g hs -n replicate2_peakcalling\
\
# Run peak calling for both replicates and the input (merged)\
macs2 callpeak -t replicate1_filtered_unique.bam replicate2_filtered_unique.bam -c input_filtered_unique.bam -g hs -n merged_peakcalling\
\
# Step 4: Remove blacklisted regions from all .narrowPeak files (replicates and merged)\
bedtools intersect -v -a replicate1_peakcalling_peaks.narrowPeak -b ENCFF356LFX.bed.gz > replicate1_no_blacklisted_peaks.narrowPeak\
bedtools intersect -v -a replicate2_peakcalling_peaks.narrowPeak -b ENCFF356LFX.bed.gz > replicate2_no_blacklisted_peaks.narrowPeak\
bedtools intersect -v -a merged_peakcalling_peaks.narrowPeak -b ENCFF356LFX.bed.gz > merged_no_blacklisted_peaks.narrowPeak\
\
# Count peaks in the filtered .narrowPeak files\
wc -l replicate1_no_blacklisted_peaks.narrowPeak\
wc -l replicate2_no_blacklisted_peaks.narrowPeak\
wc -l merged_no_blacklisted_peaks.narrowPeak\
\
# Step 5: Overlap analysis\
# Find overlaps between replicate1 and replicate2 peaks\
bedtools intersect -a replicate1_no_blacklisted_peaks.narrowPeak -b replicate2_no_blacklisted_peaks.narrowPeak -u | wc -l\
\
# Find the closest peaks between replicate1 and replicate2 summit regions\
bedtools closest -a replicate1_summits.bed -b replicate2_summits.bed -d > closest_output.bed\
awk '$NF >= 0 && $NF <= 100' closest_output.bed > filtered_output.bed\
\
# Find overlaps between replicate1 and merged peaks\
bedtools intersect -a replicate1_no_blacklisted_peaks.narrowPeak -b merged_no_blacklisted_peaks.narrowPeak -u | wc -l\
\
# Step 6: Comparison with ENCODE data\
# Download the golden star file from ENCODE, unzip, and sort it\
gunzip ENCFF410RJD.bed.gz\
sort -k1,1 -k2,2n ENCFF410RJD.bed > sorted_ENCODE.bed\
\
# Perform intersection and remove blacklisted regions\
bedtools intersect -a replicate1_no_blacklisted_peaks.narrowPeak -b replicate2_no_blacklisted_peaks.narrowPeak -u > intersection_rep1and2_peak.narrowPeak\
bedtools intersect -v -a intersection_rep1and2_peak.narrowPeak -b ENCFF356LFX.bed.gz > intersection_no_blacklisted_peaks.narrowPeak\
\
# Compute Jaccard index to compare peak overlap\
bedtools jaccard -a merged_no_blacklisted_peaks.narrowPeak -b sorted_ENCODE.bed\
\
# Step 7: Create summit file if highest Jaccard index is found for intersection between replicates 1 and 2\
awk 'BEGIN \{OFS="\\t"\} \{print $1, $2 + $10, $2 + $10 + 1, $4, $9\}' intersection_rep1and2_no_blacklisted_peaks.narrowPeak > intersection_1_2_no_blacklisted_peaks.summits.bed\
\
# Step 8: Chromatin States\
# Use the merged dataset to intersect with chromatin state annotations\
bedtools intersect -a merged_no_blacklisted_peaks.summit.bed -b K562_ChromHMM_15states.bed -wa -wb > chromatin_states_summit_merged_noBlacklisted.txt\
\
# Step 9: GREAT Analysis\
# Perform analysis on promoters and enhancers using the provided Python script\
# Please follow the script in the Python section for counting promoters and enhancers from the GREAT table\
\
# Step 10: Visualization 1\
# Find peaks that are only found in your data and not in ENCODE data\
bedtools intersect -a merged_no_blacklisted_peaks.narrowPeak -b sorted_ENCODE.bed -v > only_by_me.narrowPeak\
\
# Sort and display the top 10 peaks based on enrichment\
sort -k7,7nr only_by_me.narrowPeak | head -10\
\
# Step 11: Visualization 2 - Only Peaks from ENCODE\
bedtools intersect -a sorted_ENCODE.bed -b merged_no_blacklisted_peaks.narrowPeak -v > only_by_Encode.narrowPeak\
sort -k7,7nr only_by_Encode.narrowPeak | head -50\
}