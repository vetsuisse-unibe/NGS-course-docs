# Publishing Your HPC Exercises on GitHub

## Creating a GitHub Account

1. Visit [github.com](https://github.com)
2. Click "Sign up" in the top right corner
3. Fill in your details:
   - Enter your email address
   - Create a password
   - Choose a username (this will be your GitHub identity)
4. Verify your email address when prompted
5. Choose the free plan (GitHub Free)

## Setting Up Git Locally

1. Install Git if you haven't already:
   ```bash
   # For Ubuntu/Debian
   sudo apt-get install git

   # For CentOS/RHEL
   sudo yum install git
   ```

2. Configure Git with your identity:
   ```bash
   git config --global user.name "Your Name"
   git config --global user.email "your_email@example.com"
   ```

## Setting Up Authentication

Choose either SSH or HTTPS authentication:

### Option 1: SSH Authentication (Recommended)

1. Generate an SSH key:
   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   # Press Enter to accept default location
   # Enter a secure passphrase (optional)
   ```

2. Add the SSH key to your GitHub account:
   ```bash
   # Display your public key
   cat ~/.ssh/id_ed25519.pub
   ```

3. Copy the displayed key and add it to GitHub:
   - Go to GitHub → Settings → SSH and GPG keys
   - Click "New SSH key"
   - Paste your key and give it a title
   - Click "Add SSH key"

### Option 2: Personal Access Token (for HTTPS)

1. Generate a token on GitHub:
   - Go to GitHub → Settings → Developer Settings → Personal Access Tokens
   - Click "Generate new token (classic)"
   - Select needed permissions (at minimum: 'repo' access)
   - Copy the generated token (you won't see it again!)

2. Store the token securely for future use

## Creating a New Repository on GitHub

1. Click the "+" in the top right corner of GitHub
2. Select "New repository"
3. Fill in repository details:
   - Name: "hpc-exercises" (or your preferred name)
   - Description: "HPC exercises using SLURM"
   - Keep it Public
   - **Do not** initialize with README, .gitignore, or license
4. Click "Create repository"

## Pushing Your Local Repository

1. Make sure you're in your exercise directory:
   ```bash
   cd hpc-exercises
   ```

2. If you haven't initialized Git yet:
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   ```

3. Connect to GitHub:

   For SSH:
   ```bash
   git remote add origin git@github.com:YOUR-USERNAME/hpc-exercises.git
   ```

   For HTTPS:
   ```bash
   git remote add origin https://github.com/YOUR-USERNAME/hpc-exercises.git
   ```

4. Push your code:
   ```bash
   git push -u origin main  # or 'master' if that's your branch name
   ```

## Verifying Your Repository

1. Visit `https://github.com/YOUR-USERNAME/hpc-exercises`
2. Check that you see:
   - All your exercise files
   - Your commit history
   - The .gitignore file
   - All tracked scripts

## Troubleshooting Common Issues

### Problem: "Remote origin already exists"
```bash
# Remove existing remote
git remote remove origin
# Add new remote
git remote add origin git@github.com:YOUR-USERNAME/hpc-exercises.git
```

### Problem: "Failed to push some refs"
```bash
# Pull first
git pull --allow-unrelated-histories origin main
# Resolve any conflicts
git push -u origin main
```

### Problem: Authentication Failed
- Verify your SSH key is added to GitHub (for SSH)
- Check your personal access token is correct (for HTTPS)
- Ensure you're using the correct remote URL format

## Best Practices

1. Keep repository clean:
   - Use .gitignore for output files
   - Don't commit sensitive data
   - Don't commit large files

2. Good commit habits:
   - Write clear commit messages
   - Commit related changes together
   - Commit regularly

3. Documentation:
   - Consider adding a README.md
   - Document special requirements
   - Include setup instructions

## Maintaining Your Repository

After initial setup, use these commands for regular updates:
```bash
# Add new or modified files
git add .

# Commit changes
git commit -m "Description of changes"

# Push to GitHub
git push
```

---
**Note:** Always replace `YOUR-USERNAME` with your actual GitHub username in commands.

**Important:** Never commit sensitive information like passwords, private keys, or personal data to your repository.