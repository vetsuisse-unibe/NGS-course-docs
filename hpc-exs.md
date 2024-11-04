### HPC cluster exercises 

#### Login into the Bioinformatics server (login8.hpc.binf.unibe.ch)
Like yesterday use the remote login extension on Visual studio code and login into the bioinformatics server with your chosen username and password. 

Open the Terminal and start the exercises. 

#### Initial Setup

Before starting the exercises, we'll set up a Git repository to track our work:
#####Create and initialize the repository
```
git init hpc-exercises
cd hpc-exercises
```
Create a .gitignore file for HPC-specific files. Type the following at the prompt

```
code .gitignore
```
This opens an empty file called .gitignore in the editor window of VSC. Now add the following lines and save the file.
```
*.out
*.err
*.txt
slurm-*.out
```
Add and commit .gitignore
```
git add .gitignore
git commit -m "Initial commit: Add .gitignore for HPC output files"
```

##### Exercise 1 
In this exercises we will create a small bash script, run it locally, then submit it as a job to Slurm using sbatch, and compare the results. We will use the cat unix command to create the script  

Create a scripts directory and track it with Git
```
mkdir scripts 
cd scripts
```
Create the test script

```
cat >test.sh <<EOL
#!/bin/bash
hostname
date
sleep 30
date
EOL
```
Make the script executable

```
ls -l test.sh 
chmod +x test.sh 
ls -l test.sh 
```
Add and commit the script

```
git add test.sh
git commit -m "Add initial test script"
```
Run locally (for demo purposes only. Real work should be submitted to Slurm to run on computing nodes.)

``` 
./test.sh 
```
Submit to Slurm
```
sbatch -p pcourseb test.sh 
squeue | grep test 
```
- sbatch is the SLURM command to submit jobs
- -u $USER  is an environment variable in Unix/Linux systems that automatically contains the username of the currently logged-in user

Questions: 
1. What is the job ID of the submitted script.
2. Where is the output of the job ? 
3. Check the Git status - are any new files created that aren't tracked?

##### Create/modify test.sh in VSCode. 
1. open a New text file with VSC
2. add the following lines in the text file 
```
#!/bin/bash
hostname
date
sleep 20
date
```
3. Save and replace the test.sh file in testDir as shown below (replace the student61 with your login)
![Image of VSC-6](vsc-7.png)

4. Save the file and commit changes:
```
git add scripts/test.sh
git commit -m "Update test script with shorter sleep time"
```

4. Submit the jo by typing the following commands in the terminal. 

```
sbatch -p pcourseb test.sh 
squeue  -u $USER 
```

##### Exercise 2: Random Number Generation

Create a new script in VSCode and track it with Git:
``` 
 #!/bin/bash
for i in {1..1000}; do echo $RANDOM >>randomNumbers.txt; done
sort -n randomNumbers.txt
```
Track the new script
```
git add scripts/test2.sh
git commit -m "Add random number generation script"
```

Submit the job using _sbatch_

```
sbatch -p pcourseb -N 1 -n 1 --mem 100 -t 2:00:00 -o test2.out -e test2.err test2.sh
```
Questions: 
1. How many files has the job submission created ? 
2. What are the contents ? 
3. Check git status - which files are untracked? Why?

##### Exercise 3: Resource Allocation in Script

A better approach is defining resource allocation inside the shell script. This way you will not need to remember for the next time and simply re-run the analysis if required.
*Create a script with embedded SLURM parameters:*
open a text file in VSC, add the below lines and save as test3.sh

```
#!/bin/bash
#SBATCH -p pcourseb # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH -t 0-2:00 # time (D-HH:MM)
#SBATCH -o test3.out # STDOUT
#SBATCH -e test3.err # STDERR


for i in {1..100}; do echo $RANDOM >> randomIntegers.txt; done
sort -n randomIntegers.txt

```
Track the new script
```
git add scripts/test3.sh
git commit -m "Add script with embedded SLURM parameters"
```
submit the job using _sbatch_
```
sbatch test3.sh 
```
##### Exercise 4: Job Control
*Create a long-running script to practice job cancellation:*
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
Track the new script
```
git add scripts/test4.sh
git commit -m "Add long-running test script"
```
Submit and cancel the job
```
sbatch test4.sh
```
Wait for job to start, then:
```
scancel <job_id>  # or scancel -u $USER
```

##### Exercise 5: Job Monitoring
*Use sacct to analyze job performance*
The sacct command can tell you information about both running jobs and finished jobs. It communicates with SLURM's database of job information and can tell you lots of useful statistics about your jobs, such as how much memory and CPU they used. When you run sacct without any arguments, it will display a summary of all completed jobs in the system. This summary may include information such as job IDs, user names, job status, start and end times, and other job-related details.
One can use the -j/--jobs flag, where it takes the job ID as the input.

Trying running sacct without parameters or with -j flag and answer the following questions:  
1. How many Jobs completed 
2. How much of memory did you request and how much was used ? Can we use this to reduce the amount of memory requested next time ? 
3. Review the Git log - how many commits have you made ? 

```
git log --oneline
```
##### Best Practices

1. Always commit your scripts before running them
2. Use meaningful commit messages
3. Don't track output files in Git
4. Keep your scripts organized in directories
5. Document significant changes in commit messages


We will continue with more SLURM jobs and track them with git in the rest of our exercises.
