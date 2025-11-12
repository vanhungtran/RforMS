# Deployment Guide

This document explains how to deploy the R for Mass Spectrometry book to GitHub Pages.

## Overview

The deployment process:
1. Renders the Quarto book to the `docs/` directory
2. Copies all HTML files from `docs/` to the root directory
3. Copies necessary supporting files (robots.txt, sitemap.xml, search.json, site_libs, etc.)
4. Commits and pushes only HTML files and necessary deployment files (excludes CSS/SCSS)
5. GitHub Pages serves the site from the root directory

## Deployment Methods

### Method 1: Automated GitHub Actions (Recommended)

The repository includes a GitHub Actions workflow (`.github/workflows/render-and-deploy.yml`) that automatically:
- Renders the book when you push to the `master` branch
- Copies HTML files from `docs/` to root
- Commits and pushes the changes

**To use:**
1. Simply push your changes to the `master` branch
2. The workflow will automatically run and deploy

### Method 2: Manual Deployment Script (Windows)

For Windows users, use the PowerShell script:

```powershell
.\deploy.ps1
```

**What it does:**
- Configures Git user and email
- Renders the Quarto book
- Copies HTML files from `docs/` to root
- Stages only HTML files and necessary deployment files (excludes CSS/SCSS)
- Commits and pushes to GitHub

### Method 3: Manual Deployment Script (Linux/Mac)

For Linux/Mac users, use the bash script:

```bash
chmod +x deploy.sh
./deploy.sh
```

**What it does:**
- Same as the PowerShell script, but for Unix-like systems

## Files Excluded from Deployment

The following files are **NOT** pushed to GitHub:
- `*.css` files (compiled CSS)
- `*.scss` files (except `custom.scss` which is the source file)
- Source files (`.qmd`, `.R`, `.bib`, etc.)
- Cache directories (`*_cache/`, `*_files/`)

## Files Included in Deployment

The following files **ARE** pushed to GitHub:
- All `*.html` files (copied from `docs/` to root)
- `robots.txt`, `sitemap.xml`, `search.json`
- `.nojekyll` (tells GitHub Pages not to use Jekyll)
- `site_libs/` directory (but excludes CSS files within it)
- `zoom-controls.html` (if present)
- `custom.scss` (source file)
- `_quarto.yml` (configuration file)

## GitHub Pages Configuration

To enable GitHub Pages:
1. Go to your repository settings on GitHub
2. Navigate to "Pages" in the left sidebar
3. Under "Source", select "Deploy from a branch"
4. Choose branch: `master`
5. Choose folder: `/ (root)`
6. Click "Save"

Your book will be available at: `https://vanhungtran.github.io/RforMS/`

## Troubleshooting

### Script fails with "quarto: command not found"
- Make sure Quarto is installed and in your PATH
- Install Quarto from: https://quarto.org/docs/get-started/

### Git authentication errors
- Make sure you have push access to the repository
- Configure Git credentials if needed:
  ```bash
  git config --global user.name "vanhungtran"
  git config --global user.email "tranhungydhcm@gmail.com"
  ```

### CSS files are still being committed
- Check that `.gitignore` includes `*.css` and `!custom.scss`
- Run `git rm --cached *.css` to remove CSS files from tracking

## Manual Steps (if scripts don't work)

If the automated scripts fail, you can deploy manually:

```bash
# 1. Render the book
quarto render

# 2. Copy HTML files
cp docs/*.html .

# 3. Copy supporting files
cp docs/robots.txt . 2>/dev/null || true
cp docs/sitemap.xml . 2>/dev/null || true
cp docs/search.json . 2>/dev/null || true
cp docs/.nojekyll . 2>/dev/null || true
cp -r docs/site_libs . 2>/dev/null || true

# 4. Stage files (excluding CSS)
git add *.html
git add robots.txt sitemap.xml search.json .nojekyll
git add site_libs/ --force
git reset site_libs/**/*.css
git add zoom-controls.html custom.scss _quarto.yml

# 5. Commit and push
git commit -m "Deploy book: Update HTML files"
git push origin master
```

## Contact

For issues or questions, contact: tranhungydhcm@gmail.com

