#!/bin/zsh
set -euo pipefail

######## CONFIG (change if needed) ########
LAB_DIR="/Users/zijiefeng/Desktop/Guo's lab"          # your working repo (unchanged unless you run "setup")
APP_SUBMODULE_NAME="qPCR_Cleaner_Analyzer"            # submodule folder name (kept as-is)
MIRROR_DIR="$HOME/Projects/guos-lab-mirror"           # mirror repo path
MIRROR_REMOTE="git@github.com:ZijieFeng-98/guos-lab-mirror.git"  # create empty repo first

# Heavy app/tool folders to exclude from mirror (comment out to include everything)
EXCLUDE_DIRS=("ImageJ.app" "Fiji.app" "Meta class/MZmine-2.52-macOS")
###########################################

say() { echo "\n==> $1"; }

need_git_repo() {
  cd "$1"
  if ! git rev-parse --show-toplevel >/dev/null 2>&1; then
    echo "❌ Not a git repo: $1"; exit 1
  fi
}

setup_lab_strategy() {
  need_git_repo "$LAB_DIR"

  say "Installing Git LFS if available"
  if command -v git-lfs >/dev/null 2>&1; then
    git lfs install || true
  else
    echo "⚠️  git-lfs not found; you can install later. LFS patterns will still be written."
  fi

  say "Writing .gitignore"
  cat > .gitignore << 'EOF'
# OS / editor
.DS_Store
Thumbs.db
.Rproj.user/
*.swp

# R / session junk
.Rhistory
.RData
.Ruserdata

# build/outputs (optional)
results/*.html
results/*.pdf

# apps/tools (never version these)
*.app/
*.exe
*.dll
*.so
EOF

  say "Writing .gitattributes (normalize + LFS patterns)"
  cat > .gitattributes << 'EOF'
* text=auto

# Text for good diffs
*.R text
*.r text
*.Rmd text
*.md text
*.csv text
*.tsv text
*.txt text
*.json text
*.yml text
*.yaml text

# Large/binary via LFS
*.rds filter=lfs diff=lfs merge=lfs -text
*.rda filter=lfs diff=lfs merge=lfs -text
*.tif filter=lfs diff=lfs merge=lfs -text
*.tiff filter=lfs diff=lfs merge=lfs -text
*.pdf filter=lfs diff=lfs merge=lfs -text
*.zip filter=lfs diff=lfs merge=lfs -text
*.gz  filter=lfs diff=lfs merge=lfs -text
*.bz2 filter=lfs diff=lfs merge=lfs -text
*.7z  filter=lfs diff=lfs merge=lfs -text
*.mp4 filter=lfs diff=lfs merge=lfs -text
*.mov filter=lfs diff=lfs merge=lfs -text
EOF

  say "Creating data/README.md (provenance template)"
  mkdir -p data
  cat > data/README.md << 'EOF'
# Data Notes
- Raw data: not stored in Git (too large). Keep on Box/Drive/server.
- Processed data: commit here if small (CSV/TSV). Use Git LFS if >50 MB.
- Regeneration: document scripts that transform raw → processed.
EOF

  say "Pre-commit hook to block files >95MB"
  HOOK=".git/hooks/pre-commit"
  mkdir -p "$(dirname "$HOOK")"
  cat > "$HOOK" << 'EOF'
#!/usr/bin/env bash
set -e
MAX=95000000
files=$(git diff --cached --name-only)
large=""
for f in $files; do
  if [ -f "$f" ]; then
    size=$(wc -c <"$f")
    if [ "$size" -gt "$MAX" ]; then
      large="$large\n$size  $f"
    fi
  fi
done
if [ -n "$large" ]; then
  echo "❌ The following staged file(s) exceed 95MB:"
  echo -e "$large"
  echo "Use Git LFS (or keep them outside Git)."
  exit 1
fi
EOF
  chmod +x "$HOOK"

  say "Staging & committing strategy files"
  git add .gitignore .gitattributes data/README.md || true
  git commit -m "Setup: gitignore, gitattributes, data README, pre-commit size guard" || echo "ℹ️  Nothing to commit."

  echo "\n✅ Setup complete in: $LAB_DIR"
  echo "Push when ready: git push"
}

mirror_full_structure() {
  # Ensure source exists
  [ -d "$LAB_DIR" ] || { echo "❌ Source not found: $LAB_DIR"; exit 1; }

  say "Preparing mirror folder: $MIRROR_DIR"
  mkdir -p "$MIRROR_DIR"

  say "Syncing from source (preserve structure; exclude .git and heavy app bundles)"
  EXCL=(--exclude ".git/" --exclude "$APP_SUBMODULE_NAME/.git/")
  for d in "${EXCLUDE_DIRS[@]}"; do EXCL+=(--exclude "$d/"); done

  rsync -a --delete --prune-empty-dirs "${EXCL[@]}" "$LAB_DIR/" "$MIRROR_DIR/"

  say "Writing README, .gitignore, .gitattributes in mirror"
  cat > "$MIRROR_DIR/README.md" << 'EOF'
# Guo's Lab (Mirror)
This repository is a *mirror* of the original "Guo's lab" working directory.
- Structure is preserved.
- Some heavy app bundles/tools may be excluded.
- For the Shiny app, see the public repo: qPCR_Cleaner_Analyzer.
This mirror is regenerated from the source and should not be edited directly in place of the source.
EOF

  cat > "$MIRROR_DIR/.gitignore" << 'EOF'
.DS_Store
.Rproj.user/
EOF
  echo "* text=auto" > "$MIRROR_DIR/.gitattributes"

  say "Initializing/committing mirror repo"
  cd "$MIRROR_DIR"
  if [ ! -d .git ]; then
    git init
    git checkout -b main
  fi

  if command -v git-lfs >/dev/null 2>&1; then
    git lfs install
    git lfs track "*.rds" "*.rda" "*.pdf" "*.png" "*.jpg" "*.jpeg" "*.tif" "*.tiff" "*.zip" "*.gz" "*.bz2" "*.7z" "*.mp4" "*.mov"
    git add .gitattributes .gitattributes.lock 2>/dev/null || true
  else
    echo "⚠️  git-lfs not found; large binaries may hit GitHub limits."
  fi

  git add -A
  git commit -m "Mirror snapshot of Guo's lab (full structure, selected exclusions)" || echo "ℹ️  Nothing to commit."

  if ! git remote get-url origin >/dev/null 2>&1; then
    echo "👉 Make sure an EMPTY repo exists at: $MIRROR_REMOTE"
    git remote add origin "$MIRROR_REMOTE"
  fi

  say "Pushing mirror (main)"
  git push -u origin main

  echo "\n✅ Mirror updated & pushed: $MIRROR_REMOTE"
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [setup|mirror|all]
  setup  - apply git strategy to Guo's lab (gitignore, gitattributes, pre-commit, data README)
  mirror - build/update and push the full-structure mirror repo
  all    - do both (default)
Current config:
  LAB_DIR      = $LAB_DIR
  MIRROR_DIR   = $MIRROR_DIR
  MIRROR_REMOTE= $MIRROR_REMOTE
  Excluding in mirror: ${EXCLUDE_DIRS[*]}
EOF
}

# -------- entrypoint --------
CMD="${1:-all}"
case "$CMD" in
  setup)  setup_lab_strategy ;;
  mirror) mirror_full_structure ;;
  all)    setup_lab_strategy; mirror_full_structure ;;
  *)      usage; exit 1 ;;
esac
