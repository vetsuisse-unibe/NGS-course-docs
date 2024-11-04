# HPC cluster exercises 

#### Login into the Bioinformatics server (login8.hpc.binf.unibe.ch)
Like yesterday use the remote login extension on Visual studio code and login into the bioinformatics server with your chosen username and password. 

Open the Terminal and start the exercises. 

#### Initial Setup

Before starting the exercises, we'll set up a Git repository to track our work:

*Create and initialize the repository*
```shell
mkdir course
cd course 
git init
mkdir hpc-exercises
cd hpc-exercises
```
Create a .gitignore file for HPC-specific files. Type the following at the prompt

```shell
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
```shell
git add .gitignore
git commit -m "Initial commit: Add .gitignore for HPC output files"
```

#### Exercise 1 
In this exercises we will create a small bash script, run it locally, then submit it as a job to Slurm using sbatch, and compare the results. We will use the cat unix command to create the script  

*Create a scripts directory and track it with Git*
```bash
mkdir scripts 
cd scripts
```
*Create the test script*

```bash
cat >test.sh <<EOL
#!/bin/bash
hostname
date
sleep 30
date
EOL
```
*Make the script executable*

```bash
ls -l test.sh 
chmod +x test.sh 
ls -l test.sh 
```
*Add and commit the script*

```bash
git add test.sh
git commit -m "Add initial test script"
```
Run locally (for demo purposes only. Real work should be submitted to Slurm to run on computing nodes.)

```bash
./test.sh 
```
Submit to Slurm
```bash
sbatch -p pcourseb test.sh 
squeue -u $USER
```
- sbatch is the SLURM command to submit jobs
- -u $USER  is an environment variable in Unix/Linux systems that automatically contains the username of the currently logged-in user

Questions: 
1. What is the job ID of the submitted script.
2. Where is the output of the job ? 
3. Check the Git status - are any new files created that aren't tracked? (use git status)
4. What does the command git log --oneline do ? 

##### Modify test.sh in VSCode. 
1. open test.sh file with VSC
2. edit the sleep time. 
```bash
#!/bin/bash
hostname
date
sleep 20
date
```
3. Save  the test.sh file

4. Now commit the changes:
```bash
git add test.sh
git commit -m "Update test script with shorter sleep time"
```

4. Submit the job by typing the following commands in the terminal. 

```bash
sbatch -p pcourseb test.sh 
squeue  -u $USER 
```

#### Exercise 2: Random Number Generation

Create a new script test2.sh with VSCode  with the following code and track it with git.

```bash
 #!/bin/bash
for i in {1..1000}; do echo $RANDOM >>randomNumbers.txt; done
sort -n randomNumbers.txt
```
*Track the new script*
```bash
git add test2.sh
git commit -m "Add random number generation script"
```

*Submit the job using _sbatch_*

```bash
sbatch -p pcourseb -N 1 -n 1 --mem 100 -t 2:00:00 -o test2.out -e test2.err test2.sh
```
Questions: 
1. How many files has the job submission created ? 
2. What are the contents ? 
3. Check git ls-files to see files tracked - which files are untracked? Why?

#### Exercise 3: Resource Allocation in Script

A better approach is defining resource allocation inside the shell script. This way you will not need to remember for the next time and simply re-run the analysis if required.

*Create a script with embedded SLURM parameters:*

Open a text file in VSC, add the below lines and save as test3.sh

```bash
#!/bin/bash
#SBATCH -p pcourseb 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --mem 8G 
#SBATCH -t 0-2:00 
#SBATCH -o test3.out
#SBATCH -e test3.err


for i in {1..100}; do echo $RANDOM >> randomIntegers.txt; done
sort -n randomIntegers.txt

```
*Track the new script*

```bash
git add test3.sh
git commit -m "Add script with embedded SLURM parameters"
```
submit the job using _sbatch_
```bash
sbatch test3.sh 
```

```
#!/bin/bash                     # Shebang line - tells the system this is a bash script

#SBATCH -p pcourseb            # Specifies which partition/queue to run the job on
                               # In this case, it's using the 'pcourseb' partition

#SBATCH -N 1                   # Requests 1 node for this job
                               # A node is a complete computer in the cluster

#SBATCH -n 1                   # Requests 1 CPU core/task
                               # This defines how many parallel processes to run

#SBATCH --mem 8G               # Requests 8 gigabytes of RAM for the job
                               # This is the total memory allocation

#SBATCH -t 0-2:00             # Sets the time limit for the job
                               # Format is D-HH:MM (0 days, 2 hours, 0 minutes)

#SBATCH -o test3.out          # Specifies where to write standard output (stdout)
                               # Will create a file named 'test3.out'

#SBATCH -e test3.err          # Specifies where to write standard error (stderr)
                               # Will create a file named 'test3.err'
```
#### Exercise 4: Job Control
*Create a long-running script to practice job cancellation:*

The scancel command can be used to cancel a job after its submitted. Write the following code and save the file as test4.sh. Submit the following job to SLURM. Wait for the  job to start running (status R), then cancel it prematurely using the scancel command.

```bash
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
```bash
git add test4.sh
git commit -m "Add long-running test script"
```
Submit and cancel the job
```bash
sbatch test4.sh
```
Wait for job to start, then:
```bash
scancel <job_id>  # or scancel -u $USER
```

#### Exercise 5: Job Monitoring
*Use sacct to analyze job performance*

The sacct command can tell you information about both running jobs and finished jobs. It communicates with SLURM's database of job information and can tell you lots of useful statistics about your jobs, such as how much memory and CPU they used. 

When you run sacct without any arguments, it will display a summary of all completed jobs in the system. This summary may include information such as job IDs, user names, job status, start and end times, and other job-related details.

One can use the -j/--jobs flag, where it takes the job ID as the input.

Trying running sacct without parameters or with -j flag and answer the following questions:  
1. How many Jobs completed 
2. How much of memory did you request and how much was used ? Can we use this to reduce the amount of memory requested next time ? 
3. Review the Git log - how many commits have you made ? 

```shell
git log --oneline
```
##### Best Practices

1. Always commit your scripts before running them
2. Use meaningful commit messages
3. Don't track output files in Git
4. Keep your scripts organized in directories
5. Document significant changes in commit messages


We will continue with more SLURM jobs and track them with git in the rest of our exercises.
