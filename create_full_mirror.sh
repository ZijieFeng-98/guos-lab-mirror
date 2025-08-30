#!/bin/zsh
set -euo pipefail

##### CONFIG — change these if you want #####
SRC="/Users/zijiefeng/Desktop/Guo's lab"             # source folder (UNCHANGED)
DEST="$HOME/Projects/guos-lab-mirror"                # mirror folder (new repo)
NEW_REPO_SSH="git@github.com:ZijieFeng-98/guos-lab-mirror.git"  # GitHub repo to push to
# Exclude heavy app bundles you don't need in Git (adjust/remove if you want literally everything)
EXCLUDE_DIRS=("ImageJ.app" "Fiji.app" "Meta class/MZmine-2.52-macOS")
#############################################

echo "==> Preparing mirror folder: $DEST"
mkdir -p "$DEST"

echo "==> Syncing from source (preserving structure)"
# Build rsync exclude args
EXCL=()
for d in "${EXCLUDE_DIRS[@]}"; do EXCL+=(--exclude "$d/"); done

rsync -a --delete --prune-empty-dirs "${EXCL[@]}" \
  --exclude ".git/" \
  --exclude "qPCR_Cleaner_Analyzer/.git/" \
  "$SRC/" "$DEST/"

# Optional: write a README explaining this is a mirror
cat > "$DEST/README.md" << 'EOF'
# Guo's Lab (Mirror)

This repository is a *mirror* of the original "Guo's lab" working directory.
- Structure is preserved.
- App bundles and certain heavy tools may be excluded.
- For the Shiny app, see the public submodule/repo: **qPCR_Cleaner_Analyzer**.

This mirror is regenerated from the source and should not be edited directly in place of the source.
EOF

# Minimal ignores + normalize line endings
cat > "$DEST/.gitignore" << 'EOF'
.DS_Store
.Rproj.user/
EOF
echo "* text=auto" > "$DEST/.gitattributes"

echo "==> Initializing Git repo for the mirror"
cd "$DEST"
if [ ! -d .git ]; then
  git init
  git checkout -b main
fi

# If Git LFS is available, track common binary types so pushes stay small
if command -v git-lfs >/dev/null 2>&1; then
  echo "==> Git LFS detected; configuring tracking"
  git lfs install
  git lfs track "*.rds" "*.rda" "*.pdf" "*.png" "*.jpg" "*.jpeg" "*.tif" "*.tiff" "*.zip" "*.gz" "*.bz2" "*.7z" "*.mp4" "*.mov"
  git add .gitattributes .gitattributes.lock 2>/dev/null || true
else
  echo "⚠️  Git LFS not found. Large binaries may cause push failures. Install Git LFS, then re-run."
fi

echo "==> Committing mirror snapshot"
git add -A
git commit -m "Mirror snapshot of Guo's lab (full structure, selected exclusions)" || echo "ℹ️ Nothing to commit."

echo "==> Configuring remote & pushing"
if ! git remote get-url origin >/dev/null 2>&1; then
  echo "👉 Create an EMPTY repo on GitHub named: $(basename "$NEW_REPO_SSH" .git) (no README)"
  git remote add origin "$NEW_REPO_SSH"
fi

git push -u origin main
echo "✅ Mirror pushed to: $NEW_REPO_SSH"
