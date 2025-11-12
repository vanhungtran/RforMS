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

Write-Host "✓ Additional files copied" -ForegroundColor Green
Write-Host ""

# Step 4: Stage only HTML files and necessary deployment files
Write-Host "Step 4: Staging files for commit..." -ForegroundColor Yellow

# Add HTML files
git add *.html

# Add necessary deployment files (but not CSS/SCSS)
git add robots.txt
git add sitemap.xml
git add search.json
git add .nojekyll

# Add site_libs (but exclude CSS files)
git add site_libs/ --force
git reset site_libs/**/*.css

# Add zoom-controls.html
if (Test-Path "zoom-controls.html") {
    git add zoom-controls.html
}

# Add custom.scss (source file, not compiled CSS)
if (Test-Path "custom.scss") {
    git add custom.scss
}

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

