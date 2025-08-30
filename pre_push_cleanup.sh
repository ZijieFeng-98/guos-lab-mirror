#!/bin/zsh
set -euo pipefail
LAB_DIR="/Users/zijiefeng/Desktop/Guo's lab"
cd "$LAB_DIR"

# 1) Ensure Git LFS is installed
if ! command -v git-lfs >/dev/null 2>&1; then
  echo "❌ Git LFS not found. Install it first (see instructions above), then re-run."
  exit 1
fi
git lfs install

# 2) Track common large/binary types with LFS
git lfs track "*.rds" "*.rda" "*.pdf" "*.png" "*.jpg" "*.jpeg" "*.tif" "*.tiff" "*.zip" "*.gz" "*.bz2" "*.7z" "*.mp4" "*.mov"
git add .gitattributes
git commit -m "Track large/binary files with Git LFS" || true

# 3) Rewrite history to move those types into LFS
git lfs migrate import --everything --include="*.rds,*.rda,*.pdf,*.png,*.jpg,*.jpeg,*.tif,*.tiff,*.zip,*.gz,*.bz2,*.7z,*.mp4,*.mov"

# 4) Remove heavy app bundles from history (don't version app binaries)
if ! python3 -c "import git_filter_repo" 2>/dev/null; then
  python3 -m pip install --user git-filter-repo
fi
python3 -m git_filter_repo --force --invert-paths \
  --path "ImageJ.app" \
  --path "Meta class/MZmine-2.52-macOS" || true

# 5) Ignore those bundles going forward
printf "ImageJ.app/\nMeta class/MZmine-2.52-macOS/\n" >> .gitignore
git add .gitignore
git commit -m "Ignore app bundles" || true

echo "✅ Cleanup done. Now push rewritten history:"
echo "    git push -u origin main --force"
