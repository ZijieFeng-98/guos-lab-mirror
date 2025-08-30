#!/bin/zsh
set -euo pipefail

##### CONFIG — change these if you want #####
SRC="/Users/zijiefeng/Desktop/Guo's lab"      # your source (unchanged)
DEST="$HOME/Projects/guos-lab-R-only"         # new mirror directory
NEW_REPO_SSH="git@github.com:ZijieFeng-98/guos-lab-R-only.git"  # GitHub repo to push to
EXCLUDE_SUBMODULE="qPCR_Cleaner_Analyzer"     # we'll skip the submodule (it's already its own repo)
#############################################

# 0) Prepare destination
mkdir -p "$DEST"
rm -rf "$DEST"/*
echo "==> Exporting R-related files into: $DEST"

# 1) rsync ONLY allowed file types, preserve folders, skip everything else
#    Keep: R scripts/notebooks/projects + light data/docs/config
rsync -av --prune-empty-dirs \
  --include '*/' \
  --include '*.R' --include '*.r' --include '*.Rmd' --include '*.Rproj' \
  --include '*.csv' --include '*.tsv' --include '*.txt' --include '*.md' \
  --include '*.json' --include '*.yml' --include '*.yaml' --include '*.ini' \
  --exclude "${EXCLUDE_SUBMODULE}/" \
  --exclude '*~' \
  --exclude '*' \
  "$SRC/" "$DEST/"

# 2) Add a README stub that points to the full repo and your app
cat > "$DEST/README.md" << 'EOF'
# R-only Mirror of Guo's Lab

This repository mirrors **only** R-related files (code, notebooks, light data, and config)
from the private working directory "Guo's lab".

- Full working repo (everything, including large files): kept privately by the author.
- Shiny app: see the public repo **qPCR_Cleaner_Analyzer**.

This mirror is regenerated from the source folder and is read-only relative to the original.
EOF

# 3) Make a lightweight .gitignore and .gitattributes (text normalization)
cat > "$DEST/.gitignore" << 'EOF'
.DS_Store
.Rproj.user/
EOF

echo "* text=auto" > "$DEST/.gitattributes"

# 4) Init / commit / set remote
cd "$DEST"
if [ ! -d .git ]; then
  git init
  git checkout -b main
fi

git add -A
git commit -m "Export R-only snapshot from Guo's lab"

# 5) Set remote if missing and push
if ! git remote get-url origin >/dev/null 2>&1; then
  echo "👉 Create an EMPTY repo on GitHub named: $(basename "$NEW_REPO_SSH" .git) (no README), or ensure it exists."
  git remote add origin "$NEW_REPO_SSH"
fi

git push -u origin main
echo "✅ Done. R-only mirror pushed to: $NEW_REPO_SSH"
