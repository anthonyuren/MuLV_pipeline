#!/bin/bash 
#PBS -N trimming
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:mem=32gb
#PBS -M marian.dore@imperial.ac.uk
#PBS -m ea

#trimming (paried ends) of fq files (removing of adaptors, first is read1, second is read2)

cd /csc/analysis/Cscbioinf/mdore/March17_Uren/March17/
module load cutadapt/1.2.1
module load trim-galore/0.4.0

mkdir Trimmed

for a in */*COMBINED_R1_001.fastq.gz; do trim_galore --fastqc -a TGAAAGACCCCACCTGTAGGTTTGGCAAGCTAGC -a2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --stringency 32 --paired $a ${a%_R1*}_R3_001.fastq.gz -o Trimmed; done        #2 mismatches

