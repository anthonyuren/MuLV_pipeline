#!/bin/bash 
#PBS -N Bymouse
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -M marian.dore@imperial.ac.uk
#PBS -m ea

cd /csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/

for file in `ls *_sorted_paired.bam_mate2.txt`
do
	#file2=${file%_*_sorted_paired.bam_mate2.txt}_sorted_paired_mate2.txt
	file2=${file%_*_COMBINED_merged_sorted_paired.bam_mate2.txt}_sorted_paired_mate2.txt
	mv $file $file2
done


perl /csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/matrix-bymouse.pl

