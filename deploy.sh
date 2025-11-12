#!/bin/bash
# Deployment Script for R for Mass Spectrometry Book
# This script renders the book and deploys to GitHub Pages
# Usage: ./deploy.sh

echo "========================================"
echo "R for Mass Spectrometry Book Deployment"
echo "========================================"
echo ""

# Configure Git user and email
echo "Configuring Git user and email..."
git config user.name "vanhungtran"
git config user.email "tranhungydhcm@gmail.com"
echo "✓ Git configured"
echo ""

# Step 1: Render the Quarto book
echo "Step 1: Rendering Quarto book..."
# Try to find quarto in common locations (Windows)
QUARTO_CMD=""
if command -v quarto &> /dev/null; then
    QUARTO_CMD="quarto"
elif [ -f "/c/Program Files/Quarto/bin/quarto.cmd" ]; then
    QUARTO_CMD="/c/Program Files/Quarto/bin/quarto.cmd"
elif [ -f "/c/Program Files/Quarto/bin/quarto.exe" ]; then
    QUARTO_CMD="/c/Program Files/Quarto/bin/quarto.exe"
elif [ -f "C:/Program Files/Quarto/bin/quarto.cmd" ]; then
    QUARTO_CMD="C:/Program Files/Quarto/bin/quarto.cmd"
else
    QUARTO_CMD="quarto"
fi

$QUARTO_CMD render
if [ $? -ne 0 ]; then
    echo "✗ Error rendering book!"
    exit 1
fi
echo "✓ Book rendered successfully"
echo ""

# Step 2: Copy HTML files from docs/ to root
echo "Step 2: Copying HTML files from docs/ to root..."
copied_count=0
for file in docs/*.html; do
    if [ -f "$file" ]; then
        cp "$file" .
        copied_count=$((copied_count + 1))
        echo "  Copied: $(basename "$file")"
    fi
done
echo "✓ Copied $copied_count HTML files to root"
echo ""

# Step 3: Copy other necessary files
echo "Step 3: Copying additional files..."
for file in robots.txt sitemap.xml search.json .nojekyll; do
    if [ -f "docs/$file" ]; then
        cp "docs/$file" .
        echo "  Copied: $file"
    fi
done

# Copy site_libs directory if it exists
if [ -d "docs/site_libs" ]; then
    if [ -d "site_libs" ]; then
        rm -rf site_libs
    fi
    cp -r docs/site_libs .
    echo "  Copied: site_libs/"
fi
echo "✓ Additional files copied"
echo ""

# Step 4: Stage only HTML files and necessary deployment files
echo "Step 4: Staging files for commit..."

# Add HTML files
git add *.html

# Add necessary deployment files (but not CSS/SCSS)
git add robots.txt sitemap.xml search.json .nojekyll 2>/dev/null

# Add site_libs but exclude CSS files
git add site_libs/ 2>/dev/null
git reset site_libs/**/*.css 2>/dev/null || true

# Add zoom-controls.html
if [ -f "zoom-controls.html" ]; then
    git add zoom-controls.html
fi

# Add custom.scss (source file, not compiled CSS)
if [ -f "custom.scss" ]; then
    git add custom.scss
fi

# Add _quarto.yml if modified
git add _quarto.yml

echo "✓ Files staged"
echo ""

# Step 5: Check if there are changes to commit
echo "Step 5: Checking for changes..."
status=$(git status --short)
if [ -z "$status" ]; then
    echo "⚠ No changes to commit"
    echo ""
else
    echo "Changes to commit:"
    echo "$status"
    echo ""
    
    # Step 6: Commit changes
    echo "Step 6: Committing changes..."
    commit_message="Deploy book: Update HTML files from docs/ to root $(date '+%Y-%m-%d %H:%M:%S')"
    git commit -m "$commit_message"
    if [ $? -ne 0 ]; then
        echo "✗ Error committing changes!"
        exit 1
    fi
    echo "✓ Changes committed"
    echo ""
    
    # Step 7: Push to GitHub
    echo "Step 7: Pushing to GitHub..."
    git push origin master
    if [ $? -ne 0 ]; then
        echo "✗ Error pushing to GitHub!"
        exit 1
    fi
    echo "✓ Successfully pushed to GitHub"
    echo ""
fi

# Summary
echo "========================================"
echo "Deployment Summary"
echo "========================================"
echo "Repository: https://github.com/vanhungtran/RforMS.git"
echo "Branch: master"
echo "HTML files copied: $copied_count"
echo "GitHub Pages URL: https://vanhungtran.github.io/RforMS/"
echo ""
echo "✓ Deployment completed successfully!"
echo ""

