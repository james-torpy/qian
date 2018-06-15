#!/bin/bash

#make modules to load here

module load gi/samtools/1.2

#number of cores
numcores=8

#genome directories
genomeName="hg19_ercc"
annotationName="gencode_ercc.v19.annotation"


homeDir="/home/jamtor/"
genomeDir="$homeDir/genomes/hg19_ercc/"
genomeFile="$genomeDir/$genomeName.fa"
annotationFile="$homeDir/genomes/hg19_ercc/$annotationName.gtf"
outDir="$genomeDir/star_ref"

#log directory

logDir="$genomeDir"
mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outDir
echo $outDir
echo -e

#generate the star reference files:
star_ref_line="STAR --runMode genomeGenerate \
	--sjdbGTFfile $annotationFile --sjdbOverhang 74 --genomeDir $genomeDir \
	--genomeFastaFiles $genomeFile --runThreadN $numcores --outFileNamePrefix \
	$outDir"

echo This is the star_ref_line:
echo $star_ref_line

#submit job with name 'RSEM_count_$sample' to 15 cluster cores:
qsub -N STAR_ref_$genomeName -wd $logDir -b y -j y -R y -pe smp $numcores -V $star_ref_line
