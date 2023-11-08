#### HPC cluster exercises 

##### Exercise 1 
In this exercises we will create a small bash script, run it locally, then submit it as a job to Slurm using sbatch, and compare the results.
We are going to create a file and add contents to it using the cat command. To close the file please press CTRL and D key together. 
If you are familiar with text file editors like vi or emacs you can use that. 

```
mkdir testDir 
cd testDir 
cat >test.sh 
#!/bin/bash
hostname
date
sleep 30
date
<CTRL-D>

ls -l test.sh 
chmod +x test.sh 
ls -l test.sh 
# This is just for demo purpose. Real work should be submitted 
# to Slurm to run on computing nodes. 
./test.sh 

sbatch -p courseb test.sh 
squeue | grep test 

```

Questions: 
1. What is the job ID of the submitted script.
2. Where is the output of the job ? 

###### Now repeat the exercise with visual studio code (VSC). 
1. open a New text file with VSC
2. add the following lines in the text file 
```
#!/bin/bash
hostname
date
sleep 20
date
```
3. save and replace the test.sh file in testDir as shown below (replace the student61 with your login)
![Image of VSC-6](vsc-7.png)
4. Type the following commands in the terminal. 
```
sbatch -p courseb test.sh 
squeue | grep test
```

##### Exercise 2

open a text file in VSC, add the below lines and save the file test2.sh
``` 
#!/bin/bash
for i in {1..1000}; do echo $RANDOM >>randomNumbers.txt; done
sort -n randomNumbers.txt

```
Submit the job using _sbatch_

```
sbatch -p courseb -N 1 -n 1 --mem 100 -t 2:00:00 -o test2.out -e test2.err test2.sh
```
Questions: 
1. How many files has the job submission created ? 
2. What are the contents ? 

##### Exercise 3 

A better approach is defining resource allocation inside the shell script. This way you will not need to remember for the next time and simply re-run the analysis if required.
open a text file in VSC, add the below lines and save as test3.sh
```
#!/bin/bash
#SBATCH -p courseb # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH -t 0-2:00 # time (D-HH:MM)
#SBATCH -o test3.out # STDOUT
#SBATCH -e test3.err # STDERR


for i in {1..100}; do echo $RANDOM >> randomIntegers.txt; done
sort -n randomIntegers.txt

```
submit the job using _sbatch_

```
sbatch test3.sh 
```
##### Exercise 4
The scancel command can be used to cancel a job after its submitted. Lets go ahead and resubmit the following job. Wait for the  job to start running (status R), then cancel it prematurely using the scancel command.

```
#!/bin/bash
#SBATCH -p courseb # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 10G # memory pool for all cores
#SBATCH -t 0-2:00 # time (D-HH:MM)
#SBATCH -o test4.out # STDOUT
#SBATCH -e test5.err # STDERR

hostname
date
sleep 600
date
```

```
scancel <job_id> or scancel -u <studentid>
```
##### Exercise 5
The sacct command can tell you information about both running jobs and finished jobs. It communicates with SLURM’s database of job information and can tell you lots of useful statistics about your jobs, such as how much memory and CPU they used. When you run sacct without any arguments, it will display a summary of all completed jobs in the system. This summary may include information such as job IDs, user names, job status, start and end times, and other job-related details.
One can use the -j/--jobs flag, where it takes the job ID as the input.
Trying running sacct without parameters or with -j flag and answer the following questions:  
1. How many Jobs completed 
2. How much of memory did you request and how much was used ? Can we use this to reduce the amount of memory requested next time ? 



We will continue with more SLURM jobs in the rest of our exercises. 
