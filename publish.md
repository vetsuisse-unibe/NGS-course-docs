# Publishing Bioinformatics Exercises to GitHub

## Overview
This goal of the exercise is to publish your course exercises to GitHub, starting with your HPC exercises and preparing for future work including mapping next-generation sequences and variant calling.

## Prerequisites
* Completed HPC cluster exercises with local Git repository
* GitHub account (create one at github.com if needed). *The email id you use to create the account should be the same you used in git config command*

## Exercise 1: GitHub Setup

### Create a New GitHub Repository
1. Visit github.com and log in
2. Click the '+' icon in the top right
3. Select 'New repository'
4. Name it 'bioinformatics-exercises' or 'Sequencing-course-exercises' or choose a name repository that clearly reflects its purpose or content.(Bad:my-Project, my-first-repository etc)
5. Leave it public
6. Don't initialize with README
7. Copy the repository SSH URL

You just created your first scientific notebook on GitHub!. Feel free to store all your other code and projects here too - it's like having unlimited digital lab notebooks with automatic backup. Generally, keep one GitHub repository per analysis project.


### Configure Git on HPC Cluster (Bioinformatics server)
Use remote SSH login and connect to the IBU cluster (login8.hpc.binf.unibe.ch).

### Set Up Authentication (Important!)
To push your local repository to GitHub securely, we will use SSH authentication. SSH uses a special key to verify who you are.You generate a key pair on your computer, add the public key to GitHub. GitHub uses this to verify your identity everytime you try to add/update files to the remote  repository.  It also makes it easier to update GitHub automatically. This is helpful for scripts and large datasets. 

```shell
# Generate SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"

# Display your public key (copy this output)
cat ~/.ssh/id_ed25519.pub
```
-t ed25519: Specifies the type of key (Ed25519 is a modern, secure algorithm)

When you run this, it will:
- Ask where to save the key (default location is good)
- Ask for a passphrase (optional but recommended for security)
A passphrase is like a password but typically longer and used to encrypt your SSH private key. If you don't want to set it just type _ENTER_. If you set please make sure you remember it or store it somewhere.

Then:
1. Go to GitHub.com → Settings → SSH and GPG keys → New SSH key
2. Paste your public key and save
3. Use SSH URL when adding remote (git@github.com:username/repository.git)


## Exercise 2: Connect Local to Remote

### Link Your HPC Repository
```shell
# Navigate to your existing repository
cd hpc-exercises

# Add remote repository (replace with your URL that copy it from your github repository page )
git remote add origin git@github.com:<username>/sequencencing-exercises.git

# Verify remote was added
git remote -v
```

GT
git remote -v


The `git remote -v` command shows all the remote repositories connected to your local repository, along with their URLs. The -v stands for "verbose", showing both fetch and push URLs.

## Exercise 3: Document Your Repository

### Create README.md
Create a new file called `README.md` in your repository root or project folder (hpc-exercises) using Visual Studio Code with the following content:

```markdown
# Bioinformatics Exercises

This repository contains my work for the Bioinformatics course, including:

- HPC cluster exercises
- Quality control
- Mapping next-generation sequences
- Variant calling
- [Additional topics to be added]

## Repository Structure

- `hpc-exercises/`: Contains HPC cluster practice exercises
  - `scripts/`: Shell scripts for various HPC tasks
- [Additional directories will be added as we progress]

## Environment
- All exercises are performed on the Bioinformatics server (login8.hpc.binf.unibe.ch)
- Scripts are developed and tested using SLURM job scheduler
```

### Add Documentation to Repository
```shell
git add README.md
git commit -m "Add repository documentation"
```

## Exercise 4: Publishing Your Work
 When you created the repository on GitHub, the default branch is main (after 2020). Look at the branch name displayed on the GitHub page. But in your local repository it is master. So lets rename the local repo branch as main before pushing the contents to github remote repository. 

### Push to GitHub
```shell
# To see the default branch name
git branch 
# if it is master then change it to main 
git branch -M main 

# push your local main branch to the remote Github repository
git push -u origin main

```

### Verify Publication
1. Visit your repository URL on GitHub
2. Verify all files are present
3. Check that .gitignore is working (no output files visible)
4. Review README formatting

## Exercise 5: Ongoing Workflow

### While Working on New Exercises
```shell
# Create directory for next exercise set
cd course
mkdir mapping-exercises
cd mapping-exercises

# After creating/editing files
git add .
git commit -m "Add mapping exercise solutions"
git push origin main
```

## Best Practices

### Commit Messages
* Start with a verb (Add, Update, Fix, etc.)
* Keep first line under 50 characters
* Examples:
  * "Add variant calling scripts"
  * "Update mapping parameters for better accuracy"
  * "Fix memory allocation in SLURM script"

### Repository Organization
* Maintain separate directories for different exercise types
* Use consistent naming conventions
* Keep README.md updated
* Example structure:
  ```
  bioinformatics-exercises/
  ├── README.md
  ├── hpc-exercises/
  │   ├── scripts/
  │   └── .gitignore
  ├── mapping-exercises/
  └── variant-calling/
  ```

### Privacy and Security
* Never commit sensitive data or credentials
* Review files before committing
* Maintain an appropriate .gitignore file

## Troubleshooting Guide

### Push Rejected
If your push is rejected due to remote changes:
```shell
git pull origin main
git push origin main
```

### Authentication Issues
Generate and use SSH keys:
```shell
# Generate SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"

# Display public key to copy to GitHub
cat ~/.ssh/id_ed25519.pub
```

### Common Issues
1. **Untracked Files Appearing in Git Status**
   * Check .gitignore file
   * Use `git status` to verify
   
2. **Permission Denied**
   * Verify GitHub credentials
   * Check repository permissions

## Assessment Questions

1. What command shows configured remote repositories?
2. Why should you pull before starting new work?
3. How do you verify .gitignore is working?
4. What steps should you take if sensitive information is committed?
5. What is the recommended frequency for pushing changes?

## Additional Resources

* [GitHub Documentation](https://docs.github.com)
* [Pro Git Book](https://git-scm.com/book/en/v2)
* [GitHub Guides](https://guides.github.com)