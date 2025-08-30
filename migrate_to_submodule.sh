#!/bin/zsh
set -euo pipefail

##### CONFIG #####
LAB_DIR="/Users/zijiefeng/Desktop/Guo's lab"
LAB_REMOTE="git@github.com:ZijieFeng-98/guos-lab.git"
APP_NAME="qPCR_Cleaner_Analyzer"
APP_REMOTE="git@github.com:ZijieFeng-98/qPCR_Cleaner_Analyzer.git"
##### /CONFIG #####

say_step() { echo "\n==> $1"; }

# 0) Sanity checks
say_step "Checking lab directory exists"
if [ ! -d "$LAB_DIR" ]; then
  echo "❌ Lab directory not found: $LAB_DIR"; exit 1
fi
cd "$LAB_DIR"

# 1) Ensure lab repo is initialized on 'main'
say_step "Ensuring '$LAB_DIR' is a Git repo on 'main'"
if [ ! -d ".git" ]; then
  git init
fi
git checkout -b main 2>/dev/null || git checkout main

# Minimal housekeeping to avoid EOL noise later
if [ ! -f ".gitattributes" ]; then
  echo "* text=auto" > .gitattributes
  git add .gitattributes || true
fi
if [ ! -f ".gitignore" ]; then
  printf ".DS_Store\n.Rproj.user/\n" > .gitignore
  git add .gitignore || true
fi
git commit -m "Housekeeping: add .gitattributes/.gitignore (migration prep)" || true

# 2) If an APP folder exists, back it up so submodule can be created cleanly
BACKUP_DIR=""
if [ -d "$APP_NAME" ]; then
  TS="$(date +%Y%m%d_%H%M%S)"
  BACKUP_DIR="${APP_NAME}_backup_${TS}"
  say_step "Backing up existing '$APP_NAME' to '$BACKUP_DIR'"
  mv "$APP_NAME" "$BACKUP_DIR"
fi

# 3) Add the submodule (clones the app repo into $APP_NAME)
say_step "Adding submodule: $APP_REMOTE -> $APP_NAME"
if git config -f .gitmodules --get-regexp "submodule\..*\.path" | awk '{print $2}' | grep -qx "$APP_NAME" 2>/dev/null; then
  echo "ℹ️  Submodule '$APP_NAME' already recorded in .gitmodules, skipping add."
else
  git submodule add "$APP_REMOTE" "$APP_NAME"
fi

# 4) If we had a backup with your current app files, copy them into the submodule and push
if [ -n "${BACKUP_DIR}" ]; then
  say_step "Copying your app files from backup into the submodule working tree"
  if command -v rsync >/dev/null 2>&1; then
    rsync -av --exclude ".git" "${BACKUP_DIR}/" "${APP_NAME}/"
  else
    cp -R "${BACKUP_DIR}/." "${APP_NAME}/"
    # Remove any nested .git if it existed in the backup
    [ -d "${APP_NAME}/.git" ] && rm -rf "${APP_NAME}/.git"
  fi

  say_step "Committing & pushing files inside the APP submodule"
  pushd "$APP_NAME" >/dev/null
    git checkout -b main 2>/dev/null || git checkout main
    git add -A
    git commit -m "Import app source from Guo's lab into app repo" || echo "ℹ️  Nothing to commit in app submodule."
    git push -u origin main
  popd >/dev/null
else
  echo "ℹ️  No pre-existing app content found; using repo content from '$APP_REMOTE'."
fi

# 5) Record the submodule pointer in the lab repo and push
say_step "Recording submodule pointer in lab repo & pushing"
git add .gitmodules "$APP_NAME"
git commit -m "Add $APP_NAME as git submodule" || echo "ℹ️  Nothing to commit at lab root."
# Ensure lab remote exists
if ! git remote get-url origin >/dev/null 2>&1; then
  git remote add origin "$LAB_REMOTE"
fi
git push -u origin main

# 6) Cleanup backup if created
if [ -n "${BACKUP_DIR}" ]; then
  say_step "Cleaning up backup: $BACKUP_DIR"
  rm -rf "$BACKUP_DIR"
fi

say_step "Done!"
echo "✅ Submodule '$APP_NAME' is linked and pushed."
echo "   Lab repo: $(git remote get-url origin)"
echo "   App repo: $APP_REMOTE"
echo
echo "Fresh clone instructions (any machine):"
echo "  git clone --recurse-submodules $(git remote get-url origin)"
echo "  cd \"$(basename \"$LAB_DIR\")\""
echo "  git submodule update --init --recursive"
