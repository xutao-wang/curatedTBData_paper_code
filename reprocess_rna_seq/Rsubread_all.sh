 #!/bin/bash -l

# Set SCC project
#$ -P project_name

# Specify hard time limit for the job. 
#$ -l h_rt=4:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m a

# Give job a name
#$ -N reprocess_rna_seq

# Request eight cores
#$ -pe omp 8

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o reprocess_rna_seq.qlog

#   ask for scratch space
#$ -l scratch=100G

# Submit an array job with 10 tasks
#$ -t 1-10

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $SGE_TASK_ID"
echo "=========================================================="

# Export sratoolkit to PATH
export PATH=path_to_sratoolkit/sratoolkit.2.9.6-1-centos_linux64/bin:$PATH
export PATH=path_to_sratoolkit/.local/bin:$PATH

index=$(($SGE_TASK_ID-1))

dataDir=path_to_data_directory

# Convert txt file into array
# See get_acession_file.R for how to make geo_accession.txt
mapfile -t SRSArray < $dataDir/geo_accession.txt

# Download raw .fastq files based on SRA code
# Make sure directory $dataDir/fastq exists 
parallel-fastq-dump --sra-id  ${SRSArray[$index]} --threads 16 --outdir $dataDir/fastq --split-files --gzip # For paired-end reads
## For single-end:
## parallel-fastq-dump --sra-id  ${SRSArray[$index]} --threads 16 --outdir $dataDir/fastq --gzip
sampleName=${SRSArray[$index]}


# MAKE WORKING DIR
workingDir=${sampleName}_tmp
rm -rf $TMPDIR/$workingDir
mkdir $TMPDIR/$workingDir

# get read files
reads1=$dataDir/fastq/$sampleName*_R1_*.fastq.gz
reads2=$dataDir/fastq/$sampleName*_R2_*.fastq.gz # Comment out for single end

# TRIM the reads
java -jar path_to_Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 8 $reads1  $reads2 $TMPDIR/$workingDir/reads1.paired.fastq.gz $TMPDIR/$workingDir/reads1.unpaired.fastq.gz $TMPDIR/$workingDir/reads2.paired.fastq.gz $TMPDIR/$workingDir/reads2.unpaired.fastq.gz SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36

## For single-end:
## java -jar path_to_Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -threads 8 $reads1 $TMPDIR/$workingDir/reads1.paired.fastq.gz SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36

# Rsubread
Rsubread_script=ProcessRnaSeqFeatureCounts.R # Provided in the directory
genome=path_to_genome/hg19 # Prepared by users
genes=path_to_genes/genes.gtf # Prepared by users
threads=8

module load R/4.2.1
Rscript $Rsubread_script $genome $TMPDIR/$workingDir/reads1.paired.fastq.gz $TMPDIR/$workingDir/reads2.paired.fastq.gz $genes $dataDir/rna_seq/$sampleName $threads
## For single-end:
## Rscript $Rsubread_script $genome $TMPDIR/$workingDir/reads1.paired.fastq.gz NULL $genes $dataDir/rna_seq/$sampleName $threads



rm -rf $TMPDIR/$workingDir

