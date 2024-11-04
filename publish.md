# Publishing Bioinformatics Exercises to GitHub

## Overview
This goal of the exercise is to publish your course exercises to GitHub, starting with your HPC exercises and preparing for future work including mapping next-generation sequences and variant calling.

## Prerequisites
* Completed HPC cluster exercises with local Git repository
* GitHub account (create one at github.com if needed)

## Exercise 1: GitHub Setup

### Create a New GitHub Repository
1. Visit github.com and log in
2. Click the '+' icon in the top right
3. Select 'New repository'
4. Name it 'bioinformatics-exercises'
5. Leave it public
6. Don't initialize with README
7. Copy the repository SSH URL

### Configure Git on HPC Cluster (Bioinformatics server)
Use remote SSH login and connect to the IBU cluster.
```bash
# Set your GitHub username and email
git config --global user.name "Your GitHub Username"
git config --global user.email "your.email@example.com"
```
### Set Up Authentication (Important!)
To push your local repository to GitHub securely, you’ll need to use SSH authentication.SSH uses a special key to verify who you are. This is safer and easier than typing your password every time.

It also makes it easier to update GitHub automatically. This is helpful for scripts and large datasets. 

```bash
# Generate SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"

# Display your public key (copy this output)
cat ~/.ssh/id_ed25519.pub
```
Then:
1. Go to GitHub.com → Settings → SSH and GPG keys → New SSH key
2. Paste your public key and save
3. Use SSH URL when adding remote (git@github.com:username/repository.git)

## Exercise 2: Connect Local to Remote

### Link Your HPC Repository
```bash
# Navigate to your existing repository
cd hpc-exercises

# Add remote repository (replace with your URL)
git remote add origin https://github.com/<username>/bioinformatics-exercises.git

# Verify remote was added
git remote -v
```

## Exercise 3: Document Your Repository

### Create README.md
Create a new file called `README.md` in your repository root using Visual Studio Code with the following content:

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
```bash
git add README.md
git commit -m "Add repository documentation"
```

## Exercise 4: Publishing Your Work

### Push to GitHub
```bash
# For main branch
git push -u origin main

# If your branch is named 'master'
git push -u origin master
```

### Verify Publication
1. Visit your repository URL on GitHub
2. Verify all files are present
3. Check that .gitignore is working (no output files visible)
4. Review README formatting

## Exercise 5: Ongoing Workflow

### Before Starting New Work
```bash
# Get latest changes
git pull origin main
```

### While Working on New Exercises
```bash
# Create directory for next exercise set
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
```bash
git pull origin main
git push origin main
```

### Authentication Issues
Generate and use SSH keys:
```bash
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