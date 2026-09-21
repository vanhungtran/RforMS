# Deployment Script for R for Mass Spectrometry Book
# This script renders the book and deploys to GitHub Pages
# Usage: .\deploy.ps1

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "R for Mass Spectrometry Book Deployment" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Configure Git user and email
Write-Host "Configuring Git user and email..." -ForegroundColor Yellow
git config user.name "vanhungtran"
git config user.email "tranhungydhcm@gmail.com"
Write-Host "✓ Git configured" -ForegroundColor Green
Write-Host ""

# Step 1: Render the Quarto book
Write-Host "Step 1: Rendering Quarto book..." -ForegroundColor Yellow

# R 4.5.1 on this machine cannot load the base package; point Quarto at the
# working R 4.6.1 install explicitly so it doesn't fall back to whatever R
# is first on PATH/registry.
if (Test-Path "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe") {
    $env:QUARTO_R = "C:\Program Files\R\R-4.6.1\bin\x64"
    Write-Host "Using R at $env:QUARTO_R" -ForegroundColor Gray
}

quarto render
if ($LASTEXITCODE -ne 0) {
    Write-Host "✗ Error rendering book!" -ForegroundColor Red
    exit 1
}
Write-Host "✓ Book rendered successfully" -ForegroundColor Green
Write-Host ""

# Step 2: Copy HTML files from docs/ to root
Write-Host "Step 2: Copying HTML files from docs/ to root..." -ForegroundColor Yellow
$htmlFiles = Get-ChildItem -Path "docs" -Filter "*.html" -File
$copiedCount = 0

foreach ($file in $htmlFiles) {
    Copy-Item -Path $file.FullName -Destination "." -Force
    $copiedCount++
    Write-Host "  Copied: $($file.Name)" -ForegroundColor Gray
}

Write-Host "✓ Copied $copiedCount HTML files to root" -ForegroundColor Green
Write-Host ""

# Step 2b: Remove root-level HTML files that no longer exist under docs/ --
# e.g. left over from a renamed/renumbered chapter. Without this, deleted or
# renamed chapters stay live at their old URL forever.
Write-Host "Step 2b: Removing root-level HTML files with no docs/ counterpart..." -ForegroundColor Yellow
$removedCount = 0
Get-ChildItem -Path "." -Filter "*.html" -File | ForEach-Object {
    if ($_.Name -eq "zoom-controls.html") { return }   # not a rendered chapter
    if (-not (Test-Path "docs\$($_.Name)")) {
        Remove-Item -Path $_.FullName -Force
        $removedCount++
        Write-Host "  Removed (orphaned): $($_.Name)" -ForegroundColor Gray
    }
}
Write-Host "✓ Removed $removedCount orphaned HTML file(s)" -ForegroundColor Green
Write-Host ""

# Step 3: Copy other necessary files (robots.txt, sitemap.xml, search.json, etc.)
Write-Host "Step 3: Copying additional files..." -ForegroundColor Yellow
$additionalFiles = @("robots.txt", "sitemap.xml", "search.json", ".nojekyll")
foreach ($file in $additionalFiles) {
    if (Test-Path "docs\$file") {
        Copy-Item -Path "docs\$file" -Destination "." -Force
        Write-Host "  Copied: $file" -ForegroundColor Gray
    }
}

# Copy site_libs directory if it exists
if (Test-Path "docs\site_libs") {
    if (Test-Path "site_libs") {
        Remove-Item -Path "site_libs" -Recurse -Force
    }
    Copy-Item -Path "docs\site_libs" -Destination "." -Recurse -Force
    Write-Host "  Copied: site_libs/" -ForegroundColor Gray
}

# Copy per-chapter figure/asset directories (docs/<chapter>_files/) to root
# so root-served pages don't reference missing images.
$figureDirs = Get-ChildItem -Path "docs" -Filter "*_files" -Directory
foreach ($dir in $figureDirs) {
    $dest = $dir.Name
    if (Test-Path $dest) {
        Remove-Item -Path $dest -Recurse -Force
    }
    Copy-Item -Path $dir.FullName -Destination "." -Recurse -Force
    Write-Host "  Copied: $dest/" -ForegroundColor Gray
}

# Remove root-level *_files/ asset directories with no docs/ counterpart
# (same staleness problem as the HTML files above).
Get-ChildItem -Path "." -Filter "*_files" -Directory | ForEach-Object {
    if (-not (Test-Path "docs\$($_.Name)")) {
        Remove-Item -Path $_.FullName -Recurse -Force
        Write-Host "  Removed (orphaned): $($_.Name)/" -ForegroundColor Gray
    }
}

Write-Host "✓ Additional files copied" -ForegroundColor Green
Write-Host ""

# Step 4: Stage only HTML files and necessary deployment files
Write-Host "Step 4: Staging files for commit..." -ForegroundColor Yellow

# Add HTML files. Use `git add -A -- '*.html'` (a git pathspec) rather than
# a bare wildcard: Step 2b already deleted the orphaned files from disk, so
# anything that resolves the pattern against the *filesystem* (whether that
# resolution happens in the shell or via Get-ChildItem) has nothing left to
# match for them. `-A --` tells git itself to match the pattern against the
# index, which still knows about paths that were just deleted.
git add -A -- '*.html'

# Add necessary deployment files
git add robots.txt
git add sitemap.xml
git add search.json
git add .nojekyll

# Add site_libs/ -- this is the COMPILED theme output (cosmo + custom.scss
# baked together by `quarto render` into hashed bootstrap-<hash>.min.css).
# The Remove-Item + Copy-Item above already makes this a full resync, so
# `git add` picks up both new/changed files and removed stale ones.
git add -A -- site_libs/

# Add per-chapter figure/asset directories. Same pathspec reasoning as the
# HTML files above -- Get-ChildItem can't enumerate a directory that Step 3
# already removed, so it would never hand orphaned dirs to `git add`.
# Both the bare pattern and the /** form are passed explicitly -- a bare
# '*_files' alone left orphaned nested files (e.g. .../figure-html/*.png)
# unstaged, since it doesn't reliably recurse into matched directories.
git add -A -- '*_files' '*_files/**'

# Add zoom-controls.html
if (Test-Path "zoom-controls.html") {
    git add zoom-controls.html
}

# custom.scss itself is intentionally NOT staged here: it is book source and
# stays git-ignored/local-only per this project's tracking policy (see
# CLAUDE.md). Only its compiled output in site_libs/ above is published.

# Add _quarto.yml if modified
git add _quarto.yml

Write-Host "✓ Files staged" -ForegroundColor Green
Write-Host ""

# Step 5: Check if there are changes to commit
Write-Host "Step 5: Checking for changes..." -ForegroundColor Yellow
$status = git status --short
if ([string]::IsNullOrWhiteSpace($status)) {
    Write-Host "⚠ No changes to commit" -ForegroundColor Yellow
    Write-Host ""
} else {
    Write-Host "Changes to commit:" -ForegroundColor Cyan
    Write-Host $status
    Write-Host ""
    
    # Step 6: Commit changes
    Write-Host "Step 6: Committing changes..." -ForegroundColor Yellow
    $commitMessage = "Deploy book: Update HTML files from docs/ to root ($(Get-Date -Format 'yyyy-MM-dd HH:mm:ss'))"
    git commit -m $commitMessage
    if ($LASTEXITCODE -ne 0) {
        Write-Host "✗ Error committing changes!" -ForegroundColor Red
        exit 1
    }
    Write-Host "✓ Changes committed" -ForegroundColor Green
    Write-Host ""
    
    # Step 7: Push to GitHub
    Write-Host "Step 7: Pushing to GitHub..." -ForegroundColor Yellow
    git push origin master
    if ($LASTEXITCODE -ne 0) {
        Write-Host "✗ Error pushing to GitHub!" -ForegroundColor Red
        exit 1
    }
    Write-Host "✓ Successfully pushed to GitHub" -ForegroundColor Green
    Write-Host ""
}

# Summary
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Deployment Summary" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Repository: https://github.com/vanhungtran/RforMS.git" -ForegroundColor White
Write-Host "Branch: master" -ForegroundColor White
Write-Host "HTML files copied: $copiedCount" -ForegroundColor White
Write-Host "GitHub Pages URL: https://vanhungtran.github.io/RforMS/" -ForegroundColor White
Write-Host ""
Write-Host "✓ Deployment completed successfully!" -ForegroundColor Green
Write-Host ""

