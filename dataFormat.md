# Unlocking the Secrets of FASTQ Files

1. Login to IBU Bioinformatics sever (login8.hpc.binf.unibe.ch) using ssh protocol 
2. Create a new directory called dataPreprocess in your course directory (mkdir -p course/dataPreProcess). 
3. Change directory to dataPreprocess
4. Download the following fastq files from NCBI short read archive (http://ncbi.nlm.nih.gov/sra)

```shell
wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/001/SRR1027171/SRR1027171_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/001/SRR1027171/SRR1027171_2.fastq.gz
```
If wget didn't work please copy it from here 

```shell
cp /data/courses/pcourseb/fastq/SRR1027171_1.fastq.gz .
cp /data/courses/pcourseb/fastq/SRR1027171_2.fastq.gz .
```
-	The above files belong to the study- [GSE52194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52194) 
-	 Click the above link to take you to SRA (Sequence Read Archive) where you can find more information about these two sets.

Questions: 
1. What is the name and size of the files you downloaded ?
2. Use unix less & head command to see how the header line of fastq looks like. 

```
less SRR1027171_1.fastq.gz | head
```
Questions: 
Does it look like the example you saw in the lecture ?

#### Fastqc 

We will check the quality of the fastq files using the software program called fastqc. 

Fastqc documentation is available here 
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

We use the slurm cluster managment program to run the fastqc job. 

##### Request resources using the following srun command 

First use the _srun_ command to request for resources like CPU and RAM to run the program 

```shell
srun --partition=pcourseb --cpus-per-task 4 --time=2:00:00 --mem=1G --pty /bin/bash
```

- srun                         # SLURM command to allocate and execute resources
- --partition=pcourseb         # Use the 'pcourseb' partition/queue
- --cpus-per-task 4          # Request 4 CPU cores
- --time=2:00:00             # Request 1 hour of runtime
- --mem=4G                   # Request 4 GB of memory
- --pty /bin/bash            # Start an interactive bash shell

Questions: 
Do you understand all the arguments passed to _srun_ ? 

##### running fastqc at the shell prompt interactively 
Once you are logged into one of the compute nodes, load the software module in the following manner 

```shell
 module load FastQC/0.11.9-Java-11;
 ```

 Launch fastqc with the two fastq files in the following manner 

 ```shell
  fastqc --extract SRR1027171_1.fastq.gz SRR1027171_2.fastq.gz --threads  4
  ```
Questions: 
What is the module load command ? 


#### Job Script 
The above job can also be launched using bash script on the head node. In this case _SLURM_ looks for node with required resources and launches it on the cluster. 
Please exit the interactive session before launching the job. Meaning type exit on the node your are in and get back to the master node.

```shell
exit
```

##### Now create the following script with VS code and save as 'run_fastqc.sh' 
``` 
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=4
#SBATCH --partition=pcourseb

module add vital-it
module add UHTS/Quality_control/fastqc/0.11.5;
fastqc --extract SRR1027171_1.fastq.gz SRR1027171_2.fastq.gz --threads  4 

```
Sumbit the job to the cluster 

```shell
sbatch run_fastqc.sh 
```
Questions: 
1. Is your job running ? (Hint :use _squeue_)
It should take ~5 mins for the job to finish. 
2. How many output files have been produced by fastqc ? what are the types ?

#### zip files 

First lets create a new directory to store the zip files we will be downloading from the Bioinformatics server. Create a new folder called fastqc_html in your home directory on mac. 

```
mkdir fastqc_html 
cd fastqc_html 
```
Windows users create a folder on your local PC using the same names.
Please downlooad the zip files from the following links to fastqc_html folder. 
```
https://cloud.bioinformatics.unibe.ch/index.php/s/rPMRkpdQCT9EMkw
https://cloud.bioinformatics.unibe.ch/index.php/s/SQzpZbYqWTr2xpi
```

Now from the fastqc_html directory:
1. Unzip the zip files by double clicking on the files in the file explorer
2. Double click on the SRR1027171_1_fastqc/fastqc_report.html to open in a browser
3. Double click on the SRR1027171_2_fastqc/fastqc_report.html to open in a browser

Answer the following Questions: 
1. How many reads did the files have ? 
2. Do you find any test failed ? 
3. What amount over-representation you see in the two files ?
4. What are these over-represented sequences ?

#### Filtering adapters and low quality bases
Several tools exist that can remove the adapters and filter low quality base. Examples include trimmomatic,fastx cutadapt, sickle etc. Here we will use fastp to remove the adapters and low quality bases. The fastp manual is available here 
https://github.com/OpenGene/fastp

Create the bash script in VS code and save it in the dataPreprocess dir to run fastp with the fastq files in the following manner 

``` 
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --output=qual.out
#SBATCH --error=qual.err
#SBATCH --job-name=fastp
#SBATCH --cpus-per-task=4
#SBATCH --partition=pcourseb
module add UHTS/Quality_control/fastp/0.12.5;

fastp -w 4 -q 15 -z 5 -l 50  -i SRR1027171_1.fastq.gz -I SRR1027171_2.fastq.gz -o SRR1027171_1.clean.fq.gz -O SRR1027171_2.clean.fq.gz
```


The arguments passed to cut-adapt were based on the following:  
- Trim low-quality ends from reads before adapter removal if quality is less than 15 (-q 15)
- Discard trimmed reads that are shorter than 50 bases after trimming (-l 50)
- compression level for gzip output (-z 5)
- set number of threads to 4 (-w 4) 


The best is to run the fastp algorithm using a job script. This way you are recording all the parameters and can be easily added to your methods in manuscripts for the cause reproducibility research. 

#### FastQC again 
- Now run fastqc on the cleaned fastq files. 
- Similar to the previous done exercise download the *zip files from the links below to a local directory, unzip and view the fastqc_report.html file for the results

```shell
https://cloud.bioinformatics.unibe.ch/index.php/s/zrqeADR2KKPcz7s
https://cloud.bioinformatics.unibe.ch/index.php/s/Sor65mmETSJXWTe
```
- Record the changes you see in the cleaned and trimmed reads



