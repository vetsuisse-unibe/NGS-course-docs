#### Differential expression using RNA-seq 
For this exercise we will use the datasets from the study [GSE52194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52194)
The results of the analysis of this dataset had been published in the following paper https://www.nature.com/articles/srep01689
We will use a subset only for the analysis. 

1. Fastq files can be downloaded through Gene Expression Omnibus (GEO): GSE52194 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52194)
2. The libaray layout is paired sequencing
3. Library prep protocol did not preserve information on the transcribed strand. (non-stranded)

Background: 

3 samples each from normal breast tissue and 3 different breast cancer subtypes:
- Triple negative
- HER2 positive
- Non-triple negative (Luminal A/B)

##### Tasks

1. Assess read number and quality 
2. Map the reads to the human reference genome GRCh38 
3. Count the number of reads overlapping annotated genes 
4. Test for differential gene expression between tumor subtypes 
5. Login to binfservms01 using putty

In order reduce run times again we have subset of clean reads only for chr22 at _/data/courses/course32/RNA-seq/reads_

#### Mapping 
Create a directory called RNA_seq and mapping under it 

```
mkdir -p RNA_seq/mapping 
cd RNA_seq/mapping 
```
Copy any of the two paired-end reads samples from original folder/data/courses/course32/RNA-seq/reads to your mapping folder in the following manner.

```
cp /data/courses/course32/RNA-seq/reads/*_R*.fastq.gz  . 
```
These RNA-seq reads are human breast cancer samples so we will need to map them to human reference genome (GRCh38). Using the [HiSat2 algorithm](https://ccb.jhu.edu/software/hisat2/manual.shtml). HiSat2 is the next development of TopHat2. HiSat2 is BWT based algorithm, the difference being that it uses 2 different types of indexing system instead of one. The two different indexes are global and local indexes which are used for exon mapping reads and reads spanning one or more junctions respectively. 

Write a job script and submit the mapping job tot he cluster using sbatch.

```
#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=16G
#SBATCH --output=hisat2.out
#SBATCH --error=hisat2.err
#SBATCH --job-name=hisat2
#SBATCH --cpus-per-task=4
#SBATCH --partition=pcourse32

module add UHTS/Aligner/hisat/2.1.0;
hisat2 -x /data/courses/course32/RNA-seq/reference/Homo_sapiens.GRCh38.dna.chromosome.22  -1 /data/courses/course32/RNA-seq/reads/HER21_chr22_R1.fastq.gz -2 /data/courses/course32/RNA-seq/reads/HER21_chr22_R2.fastq.gz -S HER21.sam -p 4
```
The output from hisat2 is a sam file which needs to be converted to a bam file and sorted by chromosome co-ordinates for all downstream analysis

#### Task 
Write a job script to convert the sam to bam file using _samtools view_ and _samtools sort_ 
It is very important to check every step of your analysis. Does the output make sense? Is the quality of the results good enough to continue with the analysis? To check the mapping, have a look at the summary statistics Hisat2 wrote to the error file. 

The mapping stats for all files are available here /data/courses/course32/RNA_seq/mappingstats. Use those files and unix command line tools to answer the following questions 
- What is the highest/lowest overall alignment rate?
- What is the minimum/maximum number of reads that aligned concordantly and to a unique location of the genome? 

#### Count reads 
To count the number of reads overlapping annotated genes we will use the featureCounts algorithm. http://bioinf.wehi.edu.au/featureCounts/

```
#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --output=featureCounts.out
#SBATCH --error=featureCounts.err
#SBATCH --job-name=featureCounts
#SBATCH --cpus-per-task=8
#SBATCH --partition=pcourse32


module add UHTS/Analysis/subread/1.6.0;

featureCounts -p -C -s 0 -T 8 -Q 10 --tmpDir /data/courses/course32/  -a  /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/reference/Homo_sapiens.GRCh38.98.gtf -t exon -g gene_id  -o output.txt  /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/HER21.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/HER22.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/HER23.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/NonTNBC1.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/NonTNBC2.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/TNBC3.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/Normal1.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/Normal2.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/Normal3.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/TNBC1.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/TNBC2.coordSorted.bam /data/courses/cancergenomics/RNAseq_CancerGenomics_2019/chr22/mapped/TNBC3.coordSorted.bam
```
The parameters used for featureCount is as follows:  
- -p the read sequence is paired-end data.
- -C If specified, the chimeric fragments (those fragments that have their two ends aligned to different chromosomes) will NOT be counted. This option should be used together with -p (or isPairedEnd in Rsubread featureCounts) 
- -s Indicate if strand-specific read counting should be performed 0 (un- stranded), 1 (stranded) and 2 (reversely stranded) 
- -T Number of the threads. 
- -Q The minimum mapping quality score a read must satisfy in order to be counted. 
- --tmpDir Directory under which intermediate files are saved (later re- moved) 
- -a annotation file to be used for counting
- -t Specify the feature type 
- -g Specify the attribute type used to group features (eg. exons) into meta-features (eg. genes) 
