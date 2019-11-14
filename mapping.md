###  Map Data 
Create a folder calling variantCalling and coping the reference file 

```
mkdir variantCalling
cd variantCalling
mkdir refIdx
cd refIdx
cp /data/courses/course32/variant_Calling/chr14.fa .
```
We will use only chr11 of the dog genome as the reference just for short computing run times,that way you finish the exercises faster 

#### Index the reference 
We will use bwa mem for mapping the reads to the reference genome. For this we index the reference genome first. 
if load the bwa module is loaded by typing bwa on the commmand line if you get the help for bwa. The options for bwa are also available at http://bio-bwa.sourceforge.net/bwa.shtml

```
module add UHTS/Aligner/bwa/0.7.17;
bwa 
```
Indexing aids the aligner to find potential alignment sites on the reference faster , which saves time during alignment. Indexing the reference can be run once and stored. The indexes can be reused for all bwa alignments. A new index  has be built if you are working with a different reference genomes. 

create a bash script for indexing the genome 

```
#!/bin/bash
# Slurm options
#SBATCH --job-name="bwaIdx"
#SBATCH --chdir=.
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH -p pcourse32

module add vital-it;
module add UHTS/Aligner/bwa/0.7.17;

bwa index -a bwtsw -p chr14.fa chr14.fa
```

#### Mapping 

We will use bwa mem algorithm for mapping. This is an algorithm that is most popular mapping but does not have a publication for it !! 

##### genome reads 
Create and directory called mapping 
Copy the sample genome reads for two  Bull Terriers dogs to the mapping directory. The Bull Terriers samples are named as: 
1. BT012 
2. BT134

```
cd ../
mkdir mapping 
cd mapping 
cp /data/courses/course32/variant_Calling/*.gz .
```

Mapping involves three steps: 
**Do not execute the following steps**
1. Map the reads to the indexed reference genome 
```
bwa mem -t 8 -M -R '@RG\tID:2019111402\tPL:illumina\tPU:HHV75DSXX.4\tCN:UBern\tLB:BT134-LIB\tSM:BT134'  chr14.fa  BT134_R1.fastq.gz BT134_R2.fastq.gz >BT134.sam
```
2. Convert the sam output of the mapping to a binary bam format using _samtools view_ command
```
samtools view -@8 -h -Sb -o BT134.bam BT134.sam 
```
3. sort the bam based on chromosome co-ordinates using _samtools sort_ command 

```
samtools sort -@8 BT134.bam BT134.sorted.bam 
```

This involves a lot of reading and writing to the hard disk which is highly time consuming hence using the _piping_ power of unix we will reduce these steps to single command in the following fashion

```
#!/bin/bash
# Slurm options
#SBATCH --mail-type=fail,end
#SBATCH --job-name="mapping"
#SBATCH --chdir=.
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=8
#SBATCH -p pcourse32

module add vital-it;
module add UHTS/Aligner/bwa/0.7.17;
module add UHTS/Analysis/samtools/1.8;

bwa mem -t 8 -M -R '@RG\tID:2019111402\tPL:illumina\tPU:HHV75DSXX.4\tCN:UBern\tLB:BT134-LIB\tSM:BT134'  ../refIdx/chr14.fa  BT134_R1.fastq.gz BT134_R2.fastq.gz | samtools sort -@8 -m 5G  -o BT134.sorted.bam -
samtools index BT134.sorted.bam
```
Notice the command ends with - which means that the input to the command _samtools sort_ from _\<STDIN>_ and not from a file. 
The indexing of a bam file is done using _samtool index_. This is required for further visualing and variant calling steps. 

#### Task 

Repeat the mapping step with the second Bull Terrier sample BT012. 

#### Visualizing the mapping 

Task : 
- Login into binfservas23 using srun 

Mapping can be visualized using several tools. We will use samtools tview and IGV browser in the exercise. 

##### tview 
Samtools implements a very simple text alignment viewer called tview. It uses different colors to display mapping quality or base quality, according to users’ choice. Its easy to view the alignent over the network due its text interface for  displaying alignments. 

_samtools tview_ takes as input the bam file and reference file 

```
module add UHTS/Analysis/samtools/1.8;
samtools tview BT134.sorted.bam chr14.fa
```
tview commands:
1. left and right arrows scroll
2. press ‘?’ in the viewer for help
3. '.' to toggle between dot and nucleotide views.
4. CTRL-h and CTRL-l do “big” scrolls
5. q to quit
6. Typing g allows you to go to a specific location, in this format chromosome:location: 
- 14:5408613 (shows reference allele is **G**  and **CTT** a insertion in the sample genome.)
- 14:5731405 (shows reference allele is **T** replaced by alternate allele **G**)
- 14:5822263 

##### IGV browser 

IGV browser is another light visualizing tool for mapping. 
In order to see the mapping in an IGV browser we need to transfer the bam and its index file to your local windows file. 

Create a directory called _bamFiles_ the G:

```
pscp student51@binfservms01.unibe.ch:/home/student27/RNA-seq/*bam* G:\IGEH\_PhD_Sequencing_2018\<student51>\bamFiles
```

Type IGV on the windows search tool and open the IGV browser. 

###### Load the reference genome 
By default, IGV loads Human hg19 or the last genome used on the browser. If the genome is not dog or canFam3, we  load the reference genome  canFam3:
![IGVimage](igv.1.png)

##### Visualizing read alignments
 IGV choose File > Load from File..., select the Bull Terrier bam file, and click OK. Note that the bam and index files must be in the same directory for IGV to load these properly.
 
 In the Navigation window (as shown below) repeat and visualize the co-ordinates as we used in samtools tview 

 ![igv image 2](igv.2.png)
 More information on IGV browsers can be obtained here 
 https://software.broadinstitute.org/software/igv/UserGuide

 Question: 
 1. What is the first histogram track of IGV ? 
 2. Try loading both the bam files of the Bull Terriers. 

 We will continue with the same bam files for Variant calling and see if the variant calling algorithms find the Variants that you visualized. 


