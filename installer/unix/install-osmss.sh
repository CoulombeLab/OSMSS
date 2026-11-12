#!/usr/bin/env bash
# OSMSS installer for macOS and Linux (terminal wizard; works with the bash 3.2 that ships with macOS).
# Started by "Install OSMSS (macOS).command" or "Install OSMSS (Linux).sh" in the OSMSS folder.
#
# Steps: install folder -> R check (and, on Linux, build tools) -> SAINTexpress (optional)
#        -> The Kinase Library (optional) -> copy files, build SAINTexpress, install R packages
#        -> launcher and desktop shortcut.
# This script downloads nothing itself: R installs packages from CRAN/Bioconductor, and the user
# downloads SAINTexpress and The Kinase Library files.

set -u

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OS="$(uname -s)"
MIN_R="4.4.0"
SAINT_NAME="SAINTexpress_v3.6.3__2018-03-09"
SAINT_URL="https://sourceforge.net/projects/saint-apms/files/SAINTexpress_v3.6.3__2018-03-09.tar.gz/download"
TKL_URL="https://github.com/TheKinaseLibrary/kinase-library/tree/master/src/kinase_library/databases/matrices"
case "$OS" in
  Darwin) CRAN_URL="https://cran.r-project.org/bin/macosx/" ;;
  *)      CRAN_URL="https://cran.r-project.org/bin/linux/" ;;
esac

TARGET="$HOME/OSMSS"
RSCRIPT=""; RVERSION=""
SAINT_SRC=""           # downloaded .tar.gz or unpacked folder
TKL_ST=""; TKL_FAV=""; TKL_Y=""
TMPDIRS=""
STAMP="$(date +%Y%m%d_%H%M%S)"

# ----------------------------------------------------------------------------------------
# Output and input helpers
# ----------------------------------------------------------------------------------------
if [ -t 1 ]; then B="$(printf '\033[1m')"; G="$(printf '\033[32m')"; R="$(printf '\033[31m')"; Y="$(printf '\033[33m')"; N="$(printf '\033[0m')"
else B=""; G=""; R=""; Y=""; N=""; fi
say()   { printf '%s\n' "$*"; }
ok()    { printf '%s%s%s\n' "$G" "$*" "$N"; }
warn()  { printf '%s%s%s\n' "$Y" "$*" "$N"; }
bad()   { printf '%s%s%s\n' "$R" "$*" "$N"; }
step()  { printf '\n%s== %s ==%s\n' "$B" "$*" "$N"; }
die()   { bad "$*"; cleanup; printf '\nPress Enter to close.'; read -r _; exit 1; }
cleanup() { for d in $TMPDIRS; do rm -rf "$d"; done; }
trap cleanup EXIT

ask_yn() {  # ask_yn "question" default(y/n) -> returns 0 for yes
  local q="$1" d="${2:-y}" a hint="[Y/n]"
  [ "$d" = "n" ] && hint="[y/N]"
  while true; do
    printf '%s %s ' "$q" "$hint"; read -r a || a=""
    a="$(printf '%s' "$a" | tr '[:upper:]' '[:lower:]')"
    [ -z "$a" ] && a="$d"
    case "$a" in y|yes) return 0 ;; n|no) return 1 ;; esac
  done
}

open_url() {
  say "  $1"
  if [ "$OS" = "Darwin" ]; then open "$1" >/dev/null 2>&1 || true
  elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$1" >/dev/null 2>&1 &
  fi
}

expand_path() {  # turn ~/x into $HOME/x and strip quotes/trailing slash
  local p="$1"
  p="${p%\"}"; p="${p#\"}"; p="${p%\'}"; p="${p#\'}"
  case "$p" in "~") p="$HOME" ;; "~/"*) p="$HOME/${p#\~/}" ;; esac
  [ "$p" != "/" ] && p="${p%/}"
  printf '%s' "$p"
}

# pick MODE PROMPT: MODE is file, files or folder. Prints one path per line (nothing if cancelled).
pick() {
  local mode="$1" prompt="$2" dl="$HOME/Downloads"
  [ -d "$dl" ] || dl="$HOME"
  if [ "$OS" = "Darwin" ] && command -v osascript >/dev/null 2>&1; then
    case "$mode" in
      file)   osascript -e "POSIX path of (choose file with prompt \"$prompt\" default location POSIX file \"$dl\")" 2>/dev/null ;;
      files)  osascript -e "set fs to choose file with prompt \"$prompt\" default location POSIX file \"$dl\" with multiple selections allowed" \
                        -e 'set out to ""' -e 'repeat with f in fs' -e 'set out to out & POSIX path of f & linefeed' -e 'end repeat' \
                        -e 'return out' 2>/dev/null ;;
      folder) osascript -e "POSIX path of (choose folder with prompt \"$prompt\" default location POSIX file \"$dl\")" 2>/dev/null ;;
    esac
    return
  fi
  if [ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ] && command -v zenity >/dev/null 2>&1; then
    case "$mode" in
      file)   zenity --file-selection --title="$prompt" --filename="$dl/" 2>/dev/null ;;
      files)  zenity --file-selection --multiple --separator=$'\n' --title="$prompt" --filename="$dl/" 2>/dev/null ;;
      folder) zenity --file-selection --directory --title="$prompt" --filename="$dl/" 2>/dev/null ;;
    esac
    return
  fi
  if [ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ] && command -v kdialog >/dev/null 2>&1; then
    case "$mode" in
      file)   kdialog --title "$prompt" --getopenfilename "$dl" 2>/dev/null ;;
      files)  kdialog --title "$prompt" --multiple --separate-output --getopenfilename "$dl" 2>/dev/null ;;
      folder) kdialog --title "$prompt" --getexistingdirectory "$dl" 2>/dev/null ;;
    esac
    return
  fi
  # no file picker: type (or drag) the path into the terminal
  local p
  say "$prompt" >&2
  if [ "$mode" = "files" ]; then
    say "  Type or drag the path of each file or ZIP, one per line; press Enter on an empty line to finish." >&2
    while printf '  path: ' >&2 && read -r p && [ -n "$p" ]; do expand_path "$p"; printf '\n'; done
  else
    printf '  path: ' >&2; read -r p; [ -n "$p" ] && expand_path "$p"
  fi
}

# ----------------------------------------------------------------------------------------
# R
# ----------------------------------------------------------------------------------------
ver_ge() {  # ver_ge A B -> true when version A >= B (numeric parts)
  local IFS=.
  local -a x=($1) y=($2)
  local i a b
  for i in 0 1 2 3; do
    a="${x[$i]:-0}"; b="${y[$i]:-0}"; a="${a//[^0-9]/}"; b="${b//[^0-9]/}"
    a=$((10#${a:-0})); b=$((10#${b:-0}))
    [ "$a" -gt "$b" ] && return 0
    [ "$a" -lt "$b" ] && return 1
  done
  return 0
}

r_version() { "$1" -e 'cat(as.character(getRversion()))' 2>/dev/null | tail -n 1; }

find_r() {
  local c v best="" bestv="0"
  for c in "$(command -v Rscript 2>/dev/null)" \
           /Library/Frameworks/R.framework/Resources/bin/Rscript /opt/homebrew/bin/Rscript /usr/local/bin/Rscript \
           /usr/bin/Rscript /opt/R/*/bin/Rscript; do
    [ -n "$c" ] && [ -x "$c" ] || continue
    v="$(r_version "$c")"; [ -n "$v" ] || continue
    if ver_ge "$v" "$bestv"; then best="$c"; bestv="$v"; fi
  done
  RSCRIPT="$best"; RVERSION="$([ -n "$best" ] && printf '%s' "$bestv")"
}

check_r() {
  step "Step 2 of 6: R"
  say "OSMSS runs in R. It needs R $MIN_R or later (R 4.6.1 was used for the OSMSS paper)."
  while true; do
    find_r
    if [ -n "$RSCRIPT" ] && ver_ge "$RVERSION" "$MIN_R"; then ok "R $RVERSION found: $RSCRIPT"; return; fi
    if [ -n "$RSCRIPT" ]; then bad "R $RVERSION was found, but OSMSS needs R $MIN_R or later."; else bad "R was not found on this computer."; fi
    say "Install R from CRAN (opening the download page):"
    open_url "$CRAN_URL"
    [ "$OS" != "Darwin" ] && say "  On Ubuntu/Debian, CRAN's page explains how to get a current R (the distribution's own R may be too old)."
    say ""
    say "  1) I downloaded and installed R - check again"
    say "  2) Enter the path to Rscript myself"
    say "  3) Quit"
    printf 'Choose 1, 2 or 3: '; read -r a
    case "$a" in
      2) printf 'Path to Rscript: '; read -r p; p="$(expand_path "$p")"
         if [ -x "$p" ]; then v="$(r_version "$p")"; if [ -n "$v" ] && ver_ge "$v" "$MIN_R"; then RSCRIPT="$p"; RVERSION="$v"; ok "R $v found: $p"; return; fi; fi
         bad "That is not a working Rscript $MIN_R or later." ;;
      3) die "Setup stopped. Run it again after installing R." ;;
    esac
  done
}

# Linux: tools and libraries needed to compile R packages from source
check_linux_build_tools() {
  [ "$OS" = "Linux" ] || return 0
  step "Build tools for R packages (Linux)"
  say "On Linux, R compiles most packages from source, which needs a compiler and a few libraries."
  local missing="" cmd=""
  if command -v dpkg >/dev/null 2>&1; then
    for p in build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libpng-dev; do
      dpkg -s "$p" >/dev/null 2>&1 || missing="$missing $p"
    done
    cmd="sudo apt-get install -y$missing"
  elif command -v rpm >/dev/null 2>&1; then
    for p in gcc gcc-c++ gcc-gfortran make libcurl-devel openssl-devel libxml2-devel fontconfig-devel freetype-devel libpng-devel; do
      rpm -q "$p" >/dev/null 2>&1 || missing="$missing $p"
    done
    if command -v dnf >/dev/null 2>&1; then cmd="sudo dnf install -y$missing"; else cmd="sudo yum install -y$missing"; fi
  else
    warn "Could not check the system packages on this distribution. If an R package fails to install, install a C/C++/Fortran compiler and the development files for curl, OpenSSL, libxml2, fontconfig, freetype and libpng."
    return 0
  fi
  if [ -z "$missing" ]; then ok "The build tools and libraries are installed."; return 0; fi
  warn "Missing:$missing"
  say "They can be installed with (needs your password):"
  say "  $cmd"
  if ask_yn "Run this command now?" y; then
    $cmd || warn "The command did not finish; R packages that need these libraries may fail to install."
  else
    warn "Skipped. Some R packages may fail to install; the summary at the end will list them."
  fi
}

# ----------------------------------------------------------------------------------------
# SAINTexpress and The Kinase Library
# ----------------------------------------------------------------------------------------
choose_saint() {
  step "Step 3 of 6: affinity purification (AP-MS) with SAINTexpress (optional)"
  say "To score AP-MS experiments with SAINTexpress, download SAINTexpress 3.6.3 from its authors"
  say "(file $SAINT_NAME.tar.gz, about 10 MB). Everything else in OSMSS works without it."
  if [ "$OS" = "Darwin" ]; then
    say "On macOS there is no ready-made program, so Setup compiles it from the download (about 10 minutes;"
    say "needs Apple's Command Line Tools, which Setup offers to install). Compiling on macOS has not been"
    say "tested by the OSMSS authors; if it fails, the rest of OSMSS is not affected."
  else
    say "On Linux, Setup uses the ready-made Linux program in the download, and compiles it only if that"
    say "program does not run on this system."
  fi
  if ! ask_yn "Set up SAINTexpress now? (No = skip; see README.md, \"SAINTexpress (AP-MS)\", to add it later)" y; then
    warn "SAINTexpress skipped."; return
  fi
  say "Opening the download page:"; open_url "$SAINT_URL"
  while true; do
    printf 'When the download has finished, press Enter to select the file. '; read -r _
    local f; f="$(pick file "Select $SAINT_NAME.tar.gz (or its unpacked folder)" | head -n 1)"
    if [ -z "$f" ] && ask_yn "No file selected. Select the unpacked folder instead?" n; then
      f="$(pick folder "Select the unpacked $SAINT_NAME folder" | head -n 1)"
    fi
    f="$(expand_path "$f")"
    if [ -n "$f" ] && { [ -d "$f" ] && [ -f "$f/Makefile" ] || tar -tzf "$f" 2>/dev/null | grep -q 'Makefile$'; }; then
      SAINT_SRC="$f"; ok "SAINTexpress download found: $f"; return
    fi
    bad "That is not the SAINTexpress download (expected $SAINT_NAME.tar.gz or its unpacked folder)."
    ask_yn "Try again? (No = skip SAINTexpress)" y || { warn "SAINTexpress skipped."; return; }
  done
}

valid_tkl() {  # a real matrix file: >= 5 kB, tab-separated, not a Git LFS pointer or a web page
  local f="$1" size first
  [ -f "$f" ] || return 1
  size="$(wc -c < "$f" | tr -d ' ')"; [ "$size" -ge 5000 ] || return 1
  first="$(head -n 1 "$f")"
  printf '%s' "$first" | grep -q "$(printf '\t')" || return 1
  printf '%s' "$first" | grep -qiE '<html|<!doctype|git-lfs' && return 1
  return 0
}

unpack_zip() {  # unpack_zip file.zip -> prints the temporary folder
  local d; d="$(mktemp -d "${TMPDIR:-/tmp}/osmss.XXXXXX")"; TMPDIRS="$TMPDIRS $d"
  if command -v unzip >/dev/null 2>&1; then unzip -q -o "$1" -d "$d" >/dev/null 2>&1
  else "$RSCRIPT" -e "utils::unzip('$1', exdir = '$d')" >/dev/null 2>&1; fi
  printf '%s' "$d"
}

collect_tkl() {  # collect_tkl path... (files, ZIPs or folders) -> sets TKL_ST / TKL_FAV / TKL_Y
  local p dir name f
  for p in "$@"; do
    [ -n "$p" ] || continue
    p="$(expand_path "$p")"
    case "$(basename "$p")" in
      ser_thr_all_log2_matrices.tsv)  valid_tkl "$p" && TKL_ST="$p";  continue ;;
      st_favorability.tsv)            valid_tkl "$p" && TKL_FAV="$p"; continue ;;
      tyrosine_all_log2_matrices.tsv) valid_tkl "$p" && TKL_Y="$p";   continue ;;
    esac
    if [ -d "$p" ]; then dir="$p"
    else case "$p" in *.zip|*.ZIP) dir="$(unpack_zip "$p")" ;; *) continue ;; esac; fi
    for name in ser_thr_all_log2_matrices.tsv st_favorability.tsv tyrosine_all_log2_matrices.tsv; do
      f="$(find "$dir" -type f -name "$name" 2>/dev/null | head -n 1)"
      [ -n "$f" ] && valid_tkl "$f" || continue
      case "$name" in ser_thr_all*) TKL_ST="$f" ;; st_fav*) TKL_FAV="$f" ;; tyrosine*) TKL_Y="$f" ;; esac
    done
  done
}

choose_tkl() {
  step "Step 4 of 6: kinase prediction with The Kinase Library (optional)"
  say "The kinase prediction tools use The Kinase Library's specificity matrices. Their license"
  say "(CC BY-NC-SA 3.0, non-commercial use) does not let OSMSS include them, so download them yourself:"
  say "  ser_thr/ser_thr_all_log2_matrices.tsv, ser_thr/st_favorability.tsv, tyrosine/tyrosine_all_log2_matrices.tsv"
  say "On GitHub, open each file and use 'Download raw file', or download the whole repository as a ZIP."
  if ! ask_yn "Add The Kinase Library files now? (No = skip; see README.md to add them later)" y; then
    warn "The Kinase Library skipped."; return
  fi
  say "Opening the GitHub folder:"; open_url "$TKL_URL"
  while true; do
    printf 'When the downloads have finished, press Enter to select the files (or the ZIP). '; read -r _
    local list; list="$(pick files "Select the Kinase Library .tsv files or the downloaded ZIP")"
    if [ -z "$list" ] && ask_yn "Nothing selected. Select a folder instead?" n; then list="$(pick folder "Select the folder holding the Kinase Library files")"; fi
    local IFS=$'\n'; collect_tkl $list; unset IFS
    say "Found:   ${TKL_ST:+ser_thr_all_log2_matrices.tsv }${TKL_FAV:+st_favorability.tsv }${TKL_Y:+tyrosine_all_log2_matrices.tsv}"
    if [ -n "$TKL_ST" ] && [ -n "$TKL_Y" ]; then
      ok "The required Kinase Library files are ready."
      [ -z "$TKL_FAV" ] && warn "(st_favorability.tsv is optional but recommended.)"
      return
    fi
    local need=""
    [ -z "$TKL_ST" ] && need="ser_thr_all_log2_matrices.tsv "
    [ -z "$TKL_Y" ]  && need="${need}tyrosine_all_log2_matrices.tsv"
    bad "Still needed: $need"
    say "(Files smaller than 5 kB or that are web pages are not accepted; on GitHub use 'Download raw file'.)"
    ask_yn "Select more files? (No = skip The Kinase Library)" y || { TKL_ST=""; TKL_FAV=""; TKL_Y=""; warn "The Kinase Library skipped."; return; }
  done
}

# ----------------------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------------------
copy_app() {
  if [ "$(cd "$SRC" && pwd -P)" = "$(cd "$TARGET" && pwd -P)" ]; then say "Installing in place; no files copied."; return; fi
  say "Copying OSMSS to $TARGET ..."
  local ex="--exclude=.git --exclude=logs --exclude=cache --exclude=prerun_cache --exclude=Troubelshooting --exclude=Troubleshoot --exclude=.claude --exclude=.Rproj.user --exclude=*.zip --exclude=restart.txt --exclude=.Rhistory --exclude=.RData"
  if command -v rsync >/dev/null 2>&1; then
    # shellcheck disable=SC2086
    rsync -a $ex --exclude="PSCP and Powershell commands.txt" "$SRC/" "$TARGET/" || die "Copying the OSMSS files failed."
  else
    (cd "$SRC" && tar -cf - --exclude=.git --exclude=logs --exclude=cache --exclude=prerun_cache --exclude=Troubelshooting \
       --exclude=Troubleshoot --exclude=.claude --exclude='*.zip' --exclude=restart.txt --exclude='PSCP and Powershell commands.txt' .) \
      | (cd "$TARGET" && tar -xf -) || die "Copying the OSMSS files failed."
  fi
  ok "Files copied."
}

saint_runs() { "$1" >/dev/null 2>&1; [ $? -lt 126 ]; }   # 126/127 = cannot execute / missing libraries

install_saint() {
  [ -n "$SAINT_SRC" ] || return 0
  step "Setting up SAINTexpress"
  local home="$TARGET/SAINT/$SAINT_NAME" log="$TARGET/logs/saint_build_$STAMP.log"
  mkdir -p "$TARGET/SAINT" "$TARGET/logs"
  if [ -d "$SAINT_SRC" ]; then
    [ "$(cd "$SAINT_SRC" && pwd -P)" = "$(mkdir -p "$home" && cd "$home" && pwd -P)" ] || cp -R "$SAINT_SRC/." "$home/"
  else
    tar -xzf "$SAINT_SRC" -C "$TARGET/SAINT" || { bad "Could not unpack $SAINT_SRC"; SAINT_STATUS="failed (could not unpack the download)"; return; }
    [ -d "$home" ] || { local top; top="$(tar -tzf "$SAINT_SRC" | head -n 1 | cut -d/ -f1)"; [ -n "$top" ] && [ -d "$TARGET/SAINT/$top" ] && mv "$TARGET/SAINT/$top" "$home"; }
  fi
  mkdir -p "$home/bin"
  local bin="$home/bin/SAINTexpress-int"
  if [ "$OS" = "Linux" ] && [ -f "$home/Precompiled_binaries/Linux64/SAINTexpress-int" ]; then
    cp "$home/Precompiled_binaries/Linux64/SAINTexpress-int" "$bin"; chmod +x "$bin"
    if saint_runs "$bin"; then ok "SAINTexpress (ready-made Linux program) installed: $bin"; SAINT_STATUS="set up (ready-made Linux program)"; return; fi
    warn "The ready-made Linux program does not run on this system; compiling instead."
    rm -f "$bin"
  fi
  if [ "$OS" = "Darwin" ] && ! xcode-select -p >/dev/null 2>&1; then
    warn "Compiling needs Apple's Command Line Tools. Starting their installer..."
    xcode-select --install >/dev/null 2>&1 || true
    printf 'Finish the Command Line Tools installation in the window that opened, then press Enter. '; read -r _
    xcode-select -p >/dev/null 2>&1 || { bad "The Command Line Tools are not installed; SAINTexpress was not compiled."; SAINT_STATUS="not compiled (Command Line Tools missing)"; return; }
  fi
  for t in make gcc g++; do command -v "$t" >/dev/null 2>&1 || { bad "'$t' is not installed; SAINTexpress was not compiled."; SAINT_STATUS="not compiled ($t missing)"; return; }; done
  say "Compiling SAINTexpress (about 10 minutes). Log: $log"
  ( cd "$home" && make -j 4 bin/SAINTexpress-int ) > "$log" 2>&1 &
  local pid=$! t0; t0="$(date +%s)"
  while kill -0 "$pid" 2>/dev/null; do printf '\r  compiling... %d s ' "$(( $(date +%s) - t0 ))"; sleep 5; done
  wait "$pid"; local rc=$?; printf '\n'
  if [ $rc -eq 0 ] && [ -x "$bin" ]; then ok "SAINTexpress compiled: $bin"; SAINT_STATUS="set up (compiled)"
  else bad "Compiling SAINTexpress failed (see $log). The rest of OSMSS is not affected."; SAINT_STATUS="compiling failed (see $log)"; fi
}

install_tkl() {
  [ -n "$TKL_ST" ] && [ -n "$TKL_Y" ] || return 0
  mkdir -p "$TARGET/data/matrices/ser_thr" "$TARGET/data/matrices/tyrosine"
  cp "$TKL_ST" "$TARGET/data/matrices/ser_thr/"
  [ -n "$TKL_FAV" ] && cp "$TKL_FAV" "$TARGET/data/matrices/ser_thr/"
  cp "$TKL_Y" "$TARGET/data/matrices/tyrosine/"
  ok "The Kinase Library files copied to data/matrices."
}

install_packages() {
  step "Installing R packages"
  if [ "$OS" = "Linux" ]; then say "This usually takes 20 to 60 minutes on Linux (packages are compiled). Keep this window open."
  else say "This usually takes 10 to 30 minutes. Keep this window open."; fi
  mkdir -p "$TARGET/logs"
  PKG_LOG="$TARGET/logs/install_packages_$STAMP.log"; REPORT="$TARGET/logs/install_report_$STAMP.txt"
  local script="$TARGET/installer/install_packages.R"; [ -f "$script" ] || script="$SRC/installer/install_packages.R"
  ( cd "$TARGET" && "$RSCRIPT" "$script" "$REPORT" ) 2>&1 | tee "$PKG_LOG"
}

write_launchers() {
  local l="$TARGET/start-osmss.sh"
  cat > "$l" <<EOF
#!/usr/bin/env bash
# Starts OSMSS in your web browser. Close this window (or press Ctrl+C) to stop it.
cd "\$(dirname "\$0")"
RSCRIPT="$RSCRIPT"
[ -x "\$RSCRIPT" ] || RSCRIPT="\$(command -v Rscript)"
[ -n "\$RSCRIPT" ] || { echo "R was not found. Install R from https://cran.r-project.org and run the OSMSS installer again."; read -r _; exit 1; }
echo "Starting OSMSS. It opens in your web browser; close this window to stop OSMSS."
"\$RSCRIPT" -e "shiny::runApp('.', launch.browser = TRUE)"
EOF
  chmod +x "$l"
  if [ "$OS" = "Darwin" ]; then
    printf '#!/usr/bin/env bash\nexec "%s"\n' "$l" > "$TARGET/Start OSMSS.command"; chmod +x "$TARGET/Start OSMSS.command"
  fi
}

make_shortcut() {
  local icons="$TARGET/installer/icons"
  if [ "$OS" = "Darwin" ]; then
    local app="$HOME/Desktop/OSMSS.app"
    rm -rf "$app"; mkdir -p "$app/Contents/MacOS" "$app/Contents/Resources"
    cat > "$app/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>OSMSS</string>
  <key>CFBundleDisplayName</key><string>OSMSS</string>
  <key>CFBundleIdentifier</key><string>org.osmss.launcher</string>
  <key>CFBundleVersion</key><string>3.1</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>OSMSS</string>
  <key>CFBundleIconFile</key><string>OSMSS</string>
</dict></plist>
EOF
    printf '#!/bin/bash\n# Opens OSMSS in a Terminal window (closing it stops OSMSS).\nopen -a Terminal "%s"\n' "$TARGET/Start OSMSS.command" > "$app/Contents/MacOS/OSMSS"
    chmod +x "$app/Contents/MacOS/OSMSS"
    [ -f "$icons/OSMSS.icns" ] && cp "$icons/OSMSS.icns" "$app/Contents/Resources/OSMSS.icns"
    touch "$app"
    SHORTCUT="$app"
  else
    local desk; desk="$(xdg-user-dir DESKTOP 2>/dev/null)"; [ -n "$desk" ] && [ -d "$desk" ] || desk="$HOME/Desktop"
    local entry="[Desktop Entry]
Type=Application
Name=OSMSS
Comment=Start OSMSS (One Stop MS Shop)
Exec=\"$TARGET/start-osmss.sh\"
Path=$TARGET
Icon=$icons/OSMSS.png
Terminal=true
Categories=Science;Biology;"
    mkdir -p "$HOME/.local/share/applications"
    printf '%s\n' "$entry" > "$HOME/.local/share/applications/osmss.desktop"
    chmod +x "$HOME/.local/share/applications/osmss.desktop"
    SHORTCUT="$HOME/.local/share/applications/osmss.desktop (applications menu)"
    if [ -d "$desk" ]; then
      printf '%s\n' "$entry" > "$desk/osmss.desktop"; chmod +x "$desk/osmss.desktop"
      command -v gio >/dev/null 2>&1 && gio set "$desk/osmss.desktop" metadata::trusted true >/dev/null 2>&1
      SHORTCUT="$desk/osmss.desktop and the applications menu"
    fi
  fi
}

summary() {
  step "Summary"
  say "Installed in: $TARGET"
  say "R: $RVERSION ($RSCRIPT)"
  local pre new failed status
  pre="$(awk -F'\t' '$1=="PREINSTALLED"{print $2}' "$REPORT" 2>/dev/null)"
  new="$(awk -F'\t' '$1=="INSTALLED"{print $2}' "$REPORT" 2>/dev/null)"
  failed="$(awk -F'\t' '$1=="FAILED"{print $2}' "$REPORT" 2>/dev/null)"
  status="$(awk -F'\t' '$1=="STATUS"{print $2}' "$REPORT" 2>/dev/null)"
  if [ -z "$status" ]; then bad "The R package step did not finish. See $PKG_LOG"
  else
    say "R packages newly installed: $(printf '%s' "$new" | awk -F, 'NF{print NF}END{if(!NR)print 0}' | head -n 1)"
    if [ -n "$pre" ]; then
      say "R packages you already had (left as they were, not updated): $(printf '%s' "$pre" | awk -F, '{print NF}')"
      say "  $(printf '%s' "$pre" | sed 's/,/, /g')"
      say "  If OSMSS misbehaves, an old package version may be the reason; update packages in R with BiocManager::install()."
    fi
    if [ -n "$failed" ]; then
      bad "PROBLEM: these packages could not be installed or loaded: $(printf '%s' "$failed" | sed 's/,/, /g')"
      say "  In R, try:  BiocManager::install(c(\"$(printf '%s' "$failed" | sed 's/,/", "/g')\"))"
      say "  Full log: $PKG_LOG"
    else ok "All R packages OSMSS needs are installed."; fi
  fi
  say "SAINTexpress (AP-MS): ${SAINT_STATUS:-not set up (can be added later; see README.md)}"
  if [ -n "$TKL_ST" ] && [ -n "$TKL_Y" ]; then say "The Kinase Library: added"; else say "The Kinase Library: not added (can be added later; see README.md)"; fi
  say "Shortcut: ${SHORTCUT:-none}"
  say "To start OSMSS, use the shortcut or run $TARGET/start-osmss.sh. OSMSS opens in your web browser;"
  say "closing its terminal window stops OSMSS."
}

# For testing: "OSMSS_SOURCE_ONLY=1 source install-osmss.sh" loads the functions without running Setup.
if [ -n "${OSMSS_SOURCE_ONLY:-}" ]; then trap - EXIT; return 0 2>/dev/null || exit 0; fi

# ----------------------------------------------------------------------------------------
# Self-test (non-interactive): checks the helpers on this machine
# ----------------------------------------------------------------------------------------
if [ "${1:-}" = "--self-test" ]; then
  find_r; say "R: ${RVERSION:-none} ${RSCRIPT}"
  ver_ge 4.6.1 4.4.0 && say "ver_ge 4.6.1>=4.4.0 ok"; ver_ge 4.3.3 4.4.0 || say "ver_ge 4.3.3<4.4.0 ok"; ver_ge 4.10.0 4.4.0 && say "ver_ge 4.10>=4.4 ok"
  shift
  [ $# -gt 0 ] && { collect_tkl "$@"; say "TKL: ST=${TKL_ST:-none} FAV=${TKL_FAV:-none} Y=${TKL_Y:-none}"; }
  exit 0
fi

# ----------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------
case "$OS" in Darwin|Linux) ;; *) die "This installer is for macOS and Linux. On Windows, use 'Install OSMSS (Windows).cmd'." ;; esac
clear 2>/dev/null || true
printf '%sOSMSS Setup%s (One Stop MS Shop, version 3.1)\n' "$B" "$N"
say "This installs OSMSS, a no-code R Shiny application for quantitative proteomics. It will:"
say "  1. copy OSMSS to a folder you choose,"
say "  2. check that R $MIN_R or later is installed,"
say "  3. optionally set up SAINTexpress, for affinity-purification (AP-MS) scoring,"
say "  4. optionally add The Kinase Library files, for kinase prediction,"
say "  5. install the R packages OSMSS needs (packages you already have are not updated),"
say "  6. create a shortcut that starts OSMSS."
say "It needs an internet connection. OSMSS is free software under the GNU GPL, version 3 or later;"
say "SAINTexpress and The Kinase Library are separate works with their own licenses."
ask_yn "Continue?" y || { say "Setup cancelled."; exit 0; }

step "Step 1 of 6: install folder"
while true; do
  printf 'Install folder [%s]: ' "$TARGET"; read -r t
  [ -n "$t" ] && TARGET="$(expand_path "$t")"
  case "$TARGET" in /*) ;; *) TARGET="$PWD/$TARGET" ;; esac
  if [ -d "$TARGET" ] && [ -n "$(ls -A "$TARGET" 2>/dev/null)" ] && [ ! -f "$TARGET/app.R" ] && [ "$(cd "$TARGET" && pwd -P)" != "$(cd "$SRC" && pwd -P)" ]; then
    ask_yn "$TARGET is not empty. Install OSMSS into it anyway?" n || continue
  fi
  if mkdir -p "$TARGET" 2>/dev/null && [ -w "$TARGET" ]; then ok "OSMSS will be installed in $TARGET"; break; fi
  bad "Cannot create or write to $TARGET. Choose another folder."
done

check_r
check_linux_build_tools
choose_saint
choose_tkl

step "Step 5 of 6: ready to install"
say "Install folder:      $TARGET"
say "R:                   $RVERSION"
say "SAINTexpress:        $([ -n "$SAINT_SRC" ] && printf 'will be set up' || printf 'skipped')"
say "The Kinase Library:  $([ -n "$TKL_ST" ] && [ -n "$TKL_Y" ] && printf 'will be added' || printf 'skipped')"
ask_yn "Install now?" y || { say "Setup cancelled."; exit 0; }

mkdir -p "$TARGET/logs"
copy_app
install_tkl
install_saint
install_packages
step "Step 6 of 6: shortcut"
write_launchers
if ask_yn "Create a desktop shortcut for OSMSS?" y; then make_shortcut; fi
summary
if ask_yn "Start OSMSS now?" n; then
  if [ "$OS" = "Darwin" ]; then open -a Terminal "$TARGET/Start OSMSS.command"; else "$TARGET/start-osmss.sh"; fi
fi
printf '\nSetup has finished. Press Enter to close.'; read -r _
