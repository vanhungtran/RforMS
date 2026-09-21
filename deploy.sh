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

# R 4.5.1 on this machine cannot load the base package; point Quarto at the
# working R 4.6.1 install explicitly so it doesn't fall back to whatever R
# is first on PATH/registry.
if [ -f "/c/Program Files/R/R-4.6.1/bin/x64/Rscript.exe" ]; then
    export QUARTO_R="/c/Program Files/R/R-4.6.1/bin/x64"
    echo "Using R at $QUARTO_R"
fi

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

# Step 2b: Remove root-level HTML files that no longer exist under docs/ --
# e.g. left over from a renamed/renumbered chapter. Without this, deleted or
# renamed chapters stay live at their old URL forever.
echo "Step 2b: Removing root-level HTML files with no docs/ counterpart..."
removed_count=0
for file in *.html; do
    [ -f "$file" ] || continue
    [ "$file" = "zoom-controls.html" ] && continue   # not a rendered chapter
    if [ ! -f "docs/$file" ]; then
        rm -f "$file"
        removed_count=$((removed_count + 1))
        echo "  Removed (orphaned): $file"
    fi
done
echo "✓ Removed $removed_count orphaned HTML file(s)"
echo ""

# Step 3: Copy other necessary files
echo "Step 3: Copying additional files..."
for file in robots.txt sitemap.xml search.json .nojekyll; do
    if [ -f "docs/$file" ]; then
        cp "docs/$file" .
        echo "  Copied: $file"
    fi
done

# Mirror site_libs/ (not just copy): this removes stale hashed CSS/JS bundles
# that no longer match the current theme/custom.scss compile, instead of
# letting them accumulate silently every time the theme changes.
if [ -d "docs/site_libs" ]; then
    powershell.exe -Command "robocopy 'docs/site_libs' 'site_libs' /MIR /NFL /NDL /NJH /NJS /R:1 /W:1; if (\$LASTEXITCODE -ge 8) { exit 1 } else { exit 0 }"
    if [ $? -ne 0 ]; then
        echo "✗ Error mirroring site_libs/!"
        exit 1
    fi
    echo "  Mirrored: site_libs/ (stale files removed, not just new ones added)"
fi

# Copy per-chapter figure/asset directories (docs/<chapter>_files/) to root
# so root-served pages don't reference missing images.
for dir in docs/*_files; do
    if [ -d "$dir" ]; then
        base=$(basename "$dir")
        rm -rf "$base"
        cp -r "$dir" .
        echo "  Copied: $base/"
    fi
done

# Remove root-level *_files/ asset directories with no docs/ counterpart
# (same staleness problem as the HTML files above).
for dir in *_files; do
    [ -d "$dir" ] || continue
    if [ ! -d "docs/$dir" ]; then
        rm -rf "$dir"
        echo "  Removed (orphaned): $dir/"
    fi
done
echo "✓ Additional files synced"
echo ""

# Step 4: Stage only HTML files and necessary deployment files
echo "Step 4: Staging files for commit..."

# Add HTML files. Use `git add -A -- '*.html'` (a git pathspec, quoted so
# the *shell* never expands it) rather than `git add ./*.html` (a shell
# glob): a shell glob only matches files still present on disk, so it can
# never stage the Step 2b deletions -- by the time it runs, the orphaned
# files are already gone and the glob silently has nothing to match.
git add -A -- '*.html'

# Add necessary deployment files
git add robots.txt sitemap.xml search.json .nojekyll 2>/dev/null

# Add site_libs/ -- this is the COMPILED theme output (cosmo + custom.scss
# baked together by `quarto render` into hashed bootstrap-<hash>.min.css).
# `-A` so Step 2b's robocopy /MIR removals (stale hash files from a past
# theme/custom.scss revision) are staged as deletions, not left dangling.
git add -A -- site_libs/ 2>/dev/null || true

# Add per-chapter figure/asset directories. Pathspec, not shell glob -- same
# already-deleted-so-the-glob-misses-it reasoning as the HTML files above.
# Both the bare pattern (the directory entry itself) and the /** form
# (its nested contents) are passed explicitly -- a bare '*_files' alone
# left orphaned nested files (e.g. .../figure-html/*.png) unstaged.
git add -A -- '*_files' '*_files/**' 2>/dev/null || true

# Add zoom-controls.html
if [ -f "zoom-controls.html" ]; then
    git add zoom-controls.html
fi

# custom.scss itself is intentionally NOT staged here: it is book source and
# stays git-ignored/local-only per this project's tracking policy (see
# CLAUDE.md). Only its compiled output in site_libs/ above is published.

# Add _quarto.yml if modified
git add _quarto.yml

# Add cover image if modified
git add cover.png docs/cover.png 2>/dev/null || true

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

