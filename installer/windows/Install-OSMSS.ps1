# OSMSS installer wizard for Windows (PowerShell 5.1+, Windows Forms).
# Started by "Install OSMSS (Windows).cmd" in the OSMSS folder.
#
# Steps: install folder -> R check -> SAINTexpress (optional) -> The Kinase Library (optional)
#        -> copy files and install R packages -> desktop shortcut.
# Nothing is downloaded by this script itself: R packages are installed by R from CRAN and
# Bioconductor, and SAINTexpress and The Kinase Library files are downloaded by the user.

param([switch]$SelfTest, [string]$UiTestReport)

$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
Add-Type -AssemblyName System.IO.Compression.FileSystem
[System.Windows.Forms.Application]::EnableVisualStyles()
# The console window is hidden, so report any unexpected failure in a message box.
trap {
  [void][System.Windows.Forms.MessageBox]::Show("OSMSS Setup stopped because of an error:`r`n`r`n$($_.Exception.Message)",
                                                'OSMSS Setup', 'OK', 'Error')
  exit 1
}

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------
$SourceRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$MinR       = [version]'4.4.0'
$Links = @{
  CRAN  = 'https://cran.r-project.org/bin/windows/base/'
  SAINT = 'https://sourceforge.net/projects/saint-apms/files/SAINTexpress_v3.6.3__2018-03-09.tar.gz/download'
  TKL   = 'https://github.com/TheKinaseLibrary/kinase-library/tree/master/src/kinase_library/databases/matrices'
}
$SaintFolder = 'SAINTexpress_v3.6.3__2018-03-09'
$TklFiles = [ordered]@{
  'ser_thr_all_log2_matrices.tsv'  = 'ser_thr'
  'st_favorability.tsv'            = 'ser_thr'
  'tyrosine_all_log2_matrices.tsv' = 'tyrosine'
}
$TklRequired = @('ser_thr_all_log2_matrices.tsv', 'tyrosine_all_log2_matrices.tsv')
# Folders and files never copied into an installation (local working files, private data)
$ExcludeDirs  = @('.git', 'logs', 'cache', 'prerun_cache', 'Troubelshooting', 'Troubleshoot', '.claude', '.Rproj.user')
$ExcludeFiles = @('*.zip', 'restart.txt', 'PSCP and Powershell commands.txt', '.Rhistory', '.RData')

$state = @{
  Target    = (Join-Path $env:USERPROFILE 'OSMSS')
  R         = $null      # object: Version, Rscript, Dir
  SaintDir  = $null      # folder holding SAINTexpress-int.exe
  Tkl       = @{}        # file name -> full path
  TempDirs  = New-Object System.Collections.ArrayList
  Proc      = $null
  OutLog    = $null
  ErrLog    = $null
  Report    = $null
  Started   = $null
  Done      = $false
  InstallLog = New-Object System.Collections.ArrayList
  Notes     = New-Object System.Collections.ArrayList
}

function Add-Note([string]$msg) { [void]$state.InstallLog.Add(('{0:HH:mm:ss}  {1}' -f (Get-Date), $msg)) }

# ------------------------------------------------------------------------------------------
# Helpers (no user interface)
# ------------------------------------------------------------------------------------------
function Find-R {
  $cands = New-Object System.Collections.ArrayList
  foreach ($hive in @('HKCU:', 'HKLM:')) {
    foreach ($key in @('SOFTWARE\R-core\R', 'SOFTWARE\R-core\R64', 'SOFTWARE\WOW6432Node\R-core\R')) {
      $p = "$hive\$key"
      if (Test-Path $p) {
        foreach ($sub in Get-ChildItem $p -ErrorAction SilentlyContinue) {
          $ip = (Get-ItemProperty $sub.PSPath -ErrorAction SilentlyContinue).InstallPath
          if ($ip) { [void]$cands.Add(@($sub.PSChildName, $ip)) }
        }
      }
    }
  }
  foreach ($base in @("$env:ProgramFiles\R", "$env:LOCALAPPDATA\Programs\R", "${env:ProgramFiles(x86)}\R")) {
    if ($base -and (Test-Path $base)) {
      foreach ($d in Get-ChildItem $base -Directory -Filter 'R-*' -ErrorAction SilentlyContinue) {
        [void]$cands.Add(@(($d.Name -replace '^R-', ''), $d.FullName))
      }
    }
  }
  $best = $null
  foreach ($c in $cands) {
    $rs = Join-Path $c[1] 'bin\Rscript.exe'
    if (-not (Test-Path $rs)) { continue }
    $v = $null
    try { $v = [version](($c[0] -replace '[^0-9.]', '')) } catch { continue }
    if (($null -eq $best) -or ($v -gt $best.Version)) {
      $best = [pscustomobject]@{ Version = $v; Rscript = $rs; Dir = $c[1] }
    }
  }
  return $best
}

function Get-RFromRscript([string]$rscript) {
  # Version of a user-selected Rscript.exe
  try {
    $out = & $rscript -e 'cat(as.character(getRversion()))' 2>$null
    $v = [version]($out | Select-Object -Last 1).Trim()
    return [pscustomobject]@{ Version = $v; Rscript = $rscript; Dir = (Split-Path (Split-Path (Split-Path $rscript))) }
  } catch { return $null }
}

function New-TempDir {
  $d = Join-Path $env:TEMP ('osmss_' + [guid]::NewGuid().ToString('N').Substring(0, 8))
  New-Item -ItemType Directory -Path $d | Out-Null
  [void]$state.TempDirs.Add($d)
  return $d
}

function Expand-Download([string]$path, [string[]]$wanted) {
  # Returns a folder to search: the folder itself, or a temporary folder holding the
  # unpacked archive. For tar archives only the members whose names match $wanted are
  # extracted (the SAINTexpress archive also holds thousands of source files).
  if (Test-Path $path -PathType Container) { return $path }
  $lower = $path.ToLower()
  if ($lower.EndsWith('.zip')) {
    $tmp = New-TempDir
    [System.IO.Compression.ZipFile]::ExtractToDirectory($path, $tmp)
    return $tmp
  }
  if ($lower.EndsWith('.tar.gz') -or $lower.EndsWith('.tgz') -or $lower.EndsWith('.tar')) {
    $tar = Join-Path $env:SystemRoot 'System32\tar.exe'
    if (-not (Test-Path $tar)) { throw 'This version of Windows has no tar.exe to unpack .tar.gz files. Unpack the file first (for example with 7-Zip) and select the unpacked folder.' }
    $members = & $tar -tf $path 2>$null
    if ($LASTEXITCODE -ne 0) { throw 'The file could not be read as a .tar.gz archive.' }
    $pick = @($members | Where-Object { $m = $_; ($wanted | Where-Object { $m -like "*$_" }).Count -gt 0 })
    if ($pick.Count -eq 0) { return $null }
    $tmp = New-TempDir
    & $tar -xf $path -C $tmp @pick 2>$null
    if ($LASTEXITCODE -ne 0) { throw 'tar could not unpack the file.' }
    return $tmp
  }
  return $null
}

function Find-Saint([string]$path) {
  # Returns the folder that holds the Windows SAINTexpress-int.exe, or $null.
  if ((Test-Path $path -PathType Leaf) -and ((Split-Path $path -Leaf) -ieq 'SAINTexpress-int.exe')) {
    return (Split-Path $path)
  }
  $dir = Expand-Download $path @('Windows64/SAINTexpress-int.exe', 'SAINTexpress-int.exe')
  if (-not $dir) { return $null }
  $hits = @(Get-ChildItem $dir -Recurse -Filter 'SAINTexpress-int.exe' -File -ErrorAction SilentlyContinue)
  if ($hits.Count -eq 0) { return $null }
  $win = @($hits | Where-Object { $_.FullName -match 'Windows64' })
  if ($win.Count -gt 0) { return $win[0].DirectoryName }
  return $hits[0].DirectoryName
}

function Test-TklFile([string]$file) {
  # A real matrix file, not a Git LFS pointer or an HTML page saved by mistake
  $fi = Get-Item $file
  if ($fi.Length -lt 5000) { return $false }
  $first = Get-Content $file -TotalCount 1
  return ($first -match "`t") -and ($first -notmatch '<html|<!DOCTYPE|git-lfs')
}

function Find-Tkl([string[]]$paths) {
  $found = @{}
  foreach ($p in $paths) {
    if ((Test-Path $p -PathType Leaf) -and $TklFiles.Contains((Split-Path $p -Leaf))) {
      $found[(Split-Path $p -Leaf)] = $p; continue
    }
    $dir = Expand-Download $p @($TklFiles.Keys)
    if (-not $dir) { continue }
    foreach ($name in $TklFiles.Keys) {
      $hit = @(Get-ChildItem $dir -Recurse -Filter $name -File -ErrorAction SilentlyContinue | Select-Object -First 1)
      if ($hit.Count -gt 0) { $found[$name] = $hit[0].FullName }
    }
  }
  $bad = @($found.Keys | Where-Object { -not (Test-TklFile $found[$_]) })
  foreach ($b in $bad) { $found.Remove($b) }
  return @{ Found = $found; Bad = $bad }
}

function Test-SamePath([string]$a, [string]$b) {
  $na = [System.IO.Path]::GetFullPath($a).TrimEnd('\')
  $nb = [System.IO.Path]::GetFullPath($b).TrimEnd('\')
  return $na -ieq $nb
}

function Copy-App([string]$target) {
  if (Test-SamePath $SourceRoot $target) { Add-Note 'Installing in place (the download folder is the install folder); no files copied.'; return }
  New-Item -ItemType Directory -Force -Path $target | Out-Null
  $rcArgs = @($SourceRoot, $target, '/E', '/R:1', '/W:1', '/NFL', '/NDL', '/NJH', '/NJS', '/NP', '/XD') + $ExcludeDirs + @('/XF') + $ExcludeFiles
  & robocopy @rcArgs | Out-Null
  if ($LASTEXITCODE -ge 8) { throw "Copying the OSMSS files failed (robocopy code $LASTEXITCODE)." }
  Add-Note "OSMSS files copied to $target"
}

function Install-Saint([string]$target) {
  if (-not $state.SaintDir) { return }
  $dest = Join-Path $target "SAINT\$SaintFolder\bin"
  New-Item -ItemType Directory -Force -Path $dest | Out-Null
  Copy-Item (Join-Path $state.SaintDir 'SAINTexpress-int.exe') $dest -Force
  Add-Note "SAINTexpress-int.exe copied to $dest"
}

function Install-Tkl([string]$target) {
  foreach ($name in $state.Tkl.Keys) {
    $dest = Join-Path $target ("data\matrices\" + $TklFiles[$name])
    New-Item -ItemType Directory -Force -Path $dest | Out-Null
    Copy-Item $state.Tkl[$name] $dest -Force
    Add-Note "$name copied to $dest"
  }
}

function Write-Launcher([string]$target) {
  $launcher = Join-Path $target 'Start OSMSS.cmd'
  $lines = @(
    '@echo off',
    'title OSMSS - keep this window open while you use OSMSS',
    'cd /d "%~dp0"',
    ('set "RSCRIPT=' + $state.R.Rscript + '"'),
    'if not exist "%RSCRIPT%" for /f "delims=" %%R in (''where Rscript 2^>nul'') do set "RSCRIPT=%%R"',
    'if not exist "%RSCRIPT%" (echo R was not found. Install R from https://cran.r-project.org and run the OSMSS installer again. & pause & exit /b 1)',
    'echo Starting OSMSS. It opens in your web browser; close this window to stop OSMSS.',
    '"%RSCRIPT%" -e "shiny::runApp(''.'', launch.browser = TRUE)"',
    'echo.',
    'echo OSMSS has stopped.',
    'pause'
  )
  $enc = [System.Text.Encoding]::GetEncoding([System.Globalization.CultureInfo]::CurrentCulture.TextInfo.OEMCodePage)
  [System.IO.File]::WriteAllLines($launcher, $lines, $enc)
  return $launcher
}

function New-DesktopShortcut([string]$target, [string]$launcher) {
  $desktop = [Environment]::GetFolderPath('Desktop')
  $ws  = New-Object -ComObject WScript.Shell
  $lnk = $ws.CreateShortcut((Join-Path $desktop 'OSMSS.lnk'))
  $lnk.TargetPath       = $launcher
  $lnk.WorkingDirectory = $target
  $lnk.Description      = 'Start OSMSS (One Stop MS Shop)'
  $icon = Join-Path $target 'installer\windows\OSMSS.ico'
  if (-not (Test-Path $icon)) { $icon = Join-Path $PSScriptRoot 'OSMSS.ico' }
  if (-not (Test-Path $icon)) { $icon = Join-Path $state.R.Dir 'bin\x64\Rgui.exe' }
  if (Test-Path $icon) { $lnk.IconLocation = "$icon,0" }
  $lnk.Save()
  return (Join-Path $desktop 'OSMSS.lnk')
}

function Read-Report([string]$file) {
  $r = @{}
  if (Test-Path $file) {
    foreach ($line in Get-Content $file) {
      $parts = $line -split "`t", 2
      if ($parts.Count -eq 2) { $r[$parts[0]] = $parts[1] }
    }
  }
  return $r
}

function Remove-TempDirs { foreach ($d in $state.TempDirs) { Remove-Item $d -Recurse -Force -ErrorAction SilentlyContinue } }

# Non-interactive check of the helpers (used for testing): prints what the wizard would find.
if ($SelfTest) {
  $r = Find-R
  if ($r) { "R found: $($r.Version) at $($r.Rscript) (needs >= ${MinR}: $($r.Version -ge $MinR))" } else { 'R not found' }
  "Source folder: $SourceRoot"
  return
}

# ------------------------------------------------------------------------------------------
# User interface
# ------------------------------------------------------------------------------------------
$font     = New-Object System.Drawing.Font('Segoe UI', 9.5)
$fontBold = New-Object System.Drawing.Font('Segoe UI', 13, [System.Drawing.FontStyle]::Bold)
$green = [System.Drawing.Color]::FromArgb(0, 120, 0)
$red   = [System.Drawing.Color]::FromArgb(180, 0, 0)

$form = New-Object System.Windows.Forms.Form
$form.Text = 'OSMSS Setup'
$form.Size = New-Object System.Drawing.Size(720, 540)
$form.StartPosition = 'CenterScreen'
$form.FormBorderStyle = 'FixedDialog'
$form.MaximizeBox = $false
$form.Font = $font

$header = New-Object System.Windows.Forms.Label
$header.Font = $fontBold
$header.Location = New-Object System.Drawing.Point(20, 15)
$header.Size = New-Object System.Drawing.Size(670, 30)
$form.Controls.Add($header)

$line = New-Object System.Windows.Forms.Label
$line.BorderStyle = 'Fixed3D'
$line.Location = New-Object System.Drawing.Point(0, 50)
$line.Size = New-Object System.Drawing.Size(720, 2)
$form.Controls.Add($line)

$btnBack   = New-Object System.Windows.Forms.Button
$btnNext   = New-Object System.Windows.Forms.Button
$btnCancel = New-Object System.Windows.Forms.Button
$i = 0
foreach ($b in @($btnBack, $btnNext, $btnCancel)) {
  $b.Size = New-Object System.Drawing.Size(100, 30)
  $b.Location = New-Object System.Drawing.Point((370 + 110 * $i), 455)
  $form.Controls.Add($b); $i++
}
$btnBack.Text = '< Back'; $btnNext.Text = 'Next >'; $btnCancel.Text = 'Cancel'

$pages = New-Object System.Collections.ArrayList
function New-Page([string]$title) {
  $p = New-Object System.Windows.Forms.Panel
  $p.Location = New-Object System.Drawing.Point(20, 60)
  $p.Size = New-Object System.Drawing.Size(670, 385)
  $p.Visible = $false
  $p.Tag = $title
  $form.Controls.Add($p)
  [void]$pages.Add($p)
  return $p
}
function Add-Label($page, [string]$text, [int]$y, [int]$h = 40) {
  $l = New-Object System.Windows.Forms.Label
  $l.Text = $text
  $l.Location = New-Object System.Drawing.Point(0, $y)
  $l.Size = New-Object System.Drawing.Size(665, $h)
  $page.Controls.Add($l)
  return $l
}
function Add-Link($page, [string]$text, [string]$url, [int]$y) {
  $l = New-Object System.Windows.Forms.LinkLabel
  $l.Text = $text
  $l.Location = New-Object System.Drawing.Point(0, $y)
  $l.Size = New-Object System.Drawing.Size(665, 22)
  $l.Tag = $url
  $l.Add_LinkClicked({ param($s, $e) Start-Process $s.Tag })
  $page.Controls.Add($l)
  return $l
}
function Add-Button($page, [string]$text, [int]$x, [int]$y, [int]$w = 230) {
  $b = New-Object System.Windows.Forms.Button
  $b.Text = $text
  $b.Location = New-Object System.Drawing.Point($x, $y)
  $b.Size = New-Object System.Drawing.Size($w, 30)
  $page.Controls.Add($b)
  return $b
}

# ---- Page 0: welcome ----
$p0 = New-Page 'Welcome to OSMSS Setup'
[void](Add-Label $p0 ("This wizard installs OSMSS (One Stop MS Shop), a no-code R Shiny application for the statistical " +
  "analysis of quantitative proteomics data.`r`n`r`nIt will:`r`n" +
  "  1.  copy OSMSS to a folder you choose,`r`n" +
  "  2.  check that R (version $MinR or later) is installed,`r`n" +
  "  3.  optionally set up SAINTexpress, for affinity-purification (AP-MS) scoring,`r`n" +
  "  4.  optionally add The Kinase Library files, for kinase prediction,`r`n" +
  "  5.  install the R packages OSMSS needs (packages you already have are not updated),`r`n" +
  "  6.  create a desktop shortcut that starts OSMSS.`r`n`r`n" +
  "Installing the R packages needs an internet connection and can take 10 to 30 minutes.`r`n`r`n" +
  "OSMSS is free software under the GNU General Public License, version 3 or later. SAINTexpress and " +
  "The Kinase Library are separate works with their own licenses.") 0 330)

# ---- Page 1: folder ----
$p1 = New-Page 'Choose the install folder'
[void](Add-Label $p1 'OSMSS will be installed in this folder. Your analyses, logs and downloaded reference data are also kept here.' 0 40)
$tbFolder = New-Object System.Windows.Forms.TextBox
$tbFolder.Location = New-Object System.Drawing.Point(0, 50)
$tbFolder.Size = New-Object System.Drawing.Size(540, 25)
$tbFolder.Text = $state.Target
$p1.Controls.Add($tbFolder)
$btnBrowse = Add-Button $p1 'Browse...' 550 48 110
$btnBrowse.Add_Click({
  $d = New-Object System.Windows.Forms.FolderBrowserDialog
  $d.Description = 'Choose the folder for OSMSS'
  $d.SelectedPath = Split-Path $tbFolder.Text -Parent
  if ($d.ShowDialog() -eq 'OK') {
    $tbFolder.Text = if ((Split-Path $d.SelectedPath -Leaf) -ieq 'OSMSS') { $d.SelectedPath } else { Join-Path $d.SelectedPath 'OSMSS' }
  }
})
[void](Add-Label $p1 ("A folder named OSMSS inside your user folder is a good choice. No administrator rights are needed.`r`n`r`n" +
  "If you choose a folder that already holds OSMSS, the program files are replaced and your data folder is kept.") 95 80)

# ---- Page 2: R ----
$p2 = New-Page 'Check for R'
[void](Add-Label $p2 "OSMSS runs in R. It needs R $MinR or later (R 4.6.1 was used for the OSMSS paper)." 0 25)
$lblR = Add-Label $p2 '' 35 60
$lnkCran = Add-Link $p2 'Download R for Windows from CRAN (opens your browser)' $Links.CRAN 105
[void](Add-Label $p2 "Run the R installer with its default settings, then come back here." 130 22)
$btnRecheck = Add-Button $p2 'I downloaded and installed R - check again' 0 165 320
$btnLocateR = Add-Button $p2 'Locate Rscript.exe myself...' 330 165 230
function Update-RStatus {
  if ($state.R -and $state.R.Version -ge $MinR) {
    $lblR.Text = "R $($state.R.Version) found:`r`n$($state.R.Rscript)"
    $lblR.ForeColor = $green
    $lnkCran.Visible = $false
  } elseif ($state.R) {
    $lblR.Text = "R $($state.R.Version) was found, but OSMSS needs R $MinR or later.`r`nInstall a current R from CRAN, then click 'check again'."
    $lblR.ForeColor = $red
    $lnkCran.Visible = $true
  } else {
    $lblR.Text = "R was not found on this computer.`r`nInstall R from CRAN, then click 'check again'."
    $lblR.ForeColor = $red
    $lnkCran.Visible = $true
  }
  $btnNext.Enabled = [bool]($state.R -and $state.R.Version -ge $MinR)
}
$btnRecheck.Add_Click({ $state.R = Find-R; Update-RStatus })
$btnLocateR.Add_Click({
  $d = New-Object System.Windows.Forms.OpenFileDialog
  $d.Title = 'Select Rscript.exe (in the bin folder of your R installation)'
  $d.Filter = 'Rscript.exe|Rscript.exe'
  if ($d.ShowDialog() -eq 'OK') { $state.R = Get-RFromRscript $d.FileName; Update-RStatus }
})

# ---- Page 3: SAINTexpress ----
$p3 = New-Page 'Affinity purification (AP-MS): SAINTexpress'
[void](Add-Label $p3 ("To score affinity-purification (AP-MS) experiments with SAINTexpress, download SAINTexpress 3.6.3 " +
  "from its authors (file SAINTexpress_v3.6.3__2018-03-09.tar.gz, about 10 MB). OSMSS uses its ready-made Windows " +
  "program, so nothing needs compiling. Everything else in OSMSS works without it.") 0 60)
[void](Add-Link $p3 'Download SAINTexpress 3.6.3 from SourceForge (opens your browser)' $Links.SAINT 65)
$rbSaintNow  = New-Object System.Windows.Forms.RadioButton
$rbSaintSkip = New-Object System.Windows.Forms.RadioButton
$rbSaintNow.Text  = 'Set up SAINTexpress now'
$rbSaintSkip.Text = 'Skip this step; I will add SAINTexpress later (see README.md, "SAINTexpress (AP-MS)")'
$rbSaintNow.Location  = New-Object System.Drawing.Point(0, 100);  $rbSaintNow.Size  = New-Object System.Drawing.Size(660, 24)
$rbSaintSkip.Location = New-Object System.Drawing.Point(0, 230);  $rbSaintSkip.Size = New-Object System.Drawing.Size(660, 24)
$rbSaintSkip.Checked = $true
$p3.Controls.Add($rbSaintNow); $p3.Controls.Add($rbSaintSkip)
$btnSaintFile   = Add-Button $p3 'I downloaded the file - select it...' 20 130 260
$btnSaintFolder = Add-Button $p3 'Select an unpacked folder...' 290 130 220
$lblSaint = Add-Label $p3 '' 170 50
function Set-Saint([string]$path) {
  $form.Cursor = 'WaitCursor'
  try {
    $dir = Find-Saint $path
    if ($dir) {
      $state.SaintDir = $dir; $rbSaintNow.Checked = $true
      $lblSaint.Text = "Found SAINTexpress-int.exe in the download. It will be placed in SAINT\$SaintFolder\bin."
      $lblSaint.ForeColor = $green
    } else {
      $state.SaintDir = $null
      $lblSaint.Text = 'SAINTexpress-int.exe (Windows) was not found in that file or folder. Select SAINTexpress_v3.6.3__2018-03-09.tar.gz or its unpacked folder.'
      $lblSaint.ForeColor = $red
    }
  } catch {
    $state.SaintDir = $null; $lblSaint.Text = "Could not read the file: $($_.Exception.Message)"; $lblSaint.ForeColor = $red
  } finally { $form.Cursor = 'Default' }
}
$btnSaintFile.Add_Click({
  $d = New-Object System.Windows.Forms.OpenFileDialog
  $d.Title = 'Select the SAINTexpress download'
  $d.Filter = 'SAINTexpress download (*.tar.gz, *.zip, *.exe)|*.tar.gz;*.tgz;*.zip;*.exe|All files (*.*)|*.*'
  $d.InitialDirectory = Join-Path $env:USERPROFILE 'Downloads'
  if ($d.ShowDialog() -eq 'OK') { Set-Saint $d.FileName }
})
$btnSaintFolder.Add_Click({
  $d = New-Object System.Windows.Forms.FolderBrowserDialog
  $d.Description = 'Select the unpacked SAINTexpress folder'
  if ($d.ShowDialog() -eq 'OK') { Set-Saint $d.SelectedPath }
})

# ---- Page 4: The Kinase Library ----
$p4 = New-Page 'Kinase prediction: The Kinase Library'
[void](Add-Label $p4 ("The kinase prediction tools score sites with The Kinase Library's specificity matrices. Their license " +
  "(CC BY-NC-SA 3.0, non-commercial use) does not let OSMSS include them, so download them yourself:`r`n" +
  "  - ser_thr\ser_thr_all_log2_matrices.tsv and ser_thr\st_favorability.tsv`r`n" +
  "  - tyrosine\tyrosine_all_log2_matrices.tsv`r`n" +
  "On GitHub, open each file and use 'Download raw file', or download the whole repository as a ZIP.") 0 105)
[void](Add-Link $p4 'Open The Kinase Library matrices folder on GitHub (opens your browser)' $Links.TKL 108)
$rbTklNow  = New-Object System.Windows.Forms.RadioButton
$rbTklSkip = New-Object System.Windows.Forms.RadioButton
$rbTklNow.Text  = 'Add The Kinase Library files now'
$rbTklSkip.Text = 'Skip this step; I will add them later (see README.md, "Reference data you download yourself")'
$rbTklNow.Location  = New-Object System.Drawing.Point(0, 140);  $rbTklNow.Size  = New-Object System.Drawing.Size(660, 24)
$rbTklSkip.Location = New-Object System.Drawing.Point(0, 290);  $rbTklSkip.Size = New-Object System.Drawing.Size(660, 24)
$rbTklSkip.Checked = $true
$p4.Controls.Add($rbTklNow); $p4.Controls.Add($rbTklSkip)
$btnTklFiles  = Add-Button $p4 'I downloaded them - select the files or ZIP...' 20 170 300
$btnTklFolder = Add-Button $p4 'Select a folder...' 330 170 160
$lblTkl = Add-Label $p4 '' 210 75
function Set-Tkl([string[]]$paths) {
  $form.Cursor = 'WaitCursor'
  try {
    $res = Find-Tkl $paths
    foreach ($k in $res.Found.Keys) { $state.Tkl[$k] = $res.Found[$k] }
    $have = @($TklFiles.Keys | Where-Object { $state.Tkl.ContainsKey($_) })
    $missing = @($TklRequired | Where-Object { -not $state.Tkl.ContainsKey($_) })
    $msg = 'Found: ' + $(if ($have.Count) { $have -join ', ' } else { 'none' })
    if ($res.Bad.Count) { $msg += "`r`nNot usable (too small or not a table; a GitHub page saved instead of the file?): " + ($res.Bad -join ', ') }
    if ($missing.Count) {
      $msg += "`r`nStill needed: " + ($missing -join ', '); $lblTkl.ForeColor = $red
    } else {
      $rbTklNow.Checked = $true; $lblTkl.ForeColor = $green
      if (-not $state.Tkl.ContainsKey('st_favorability.tsv')) { $msg += "`r`n(st_favorability.tsv is optional but recommended.)" }
    }
    $lblTkl.Text = $msg
  } catch {
    $lblTkl.Text = "Could not read the selection: $($_.Exception.Message)"; $lblTkl.ForeColor = $red
  } finally { $form.Cursor = 'Default' }
}
$btnTklFiles.Add_Click({
  $d = New-Object System.Windows.Forms.OpenFileDialog
  $d.Title = 'Select the Kinase Library files (.tsv) or the downloaded ZIP'
  $d.Filter = 'Kinase Library files (*.tsv, *.zip)|*.tsv;*.zip|All files (*.*)|*.*'
  $d.Multiselect = $true
  $d.InitialDirectory = Join-Path $env:USERPROFILE 'Downloads'
  if ($d.ShowDialog() -eq 'OK') { Set-Tkl $d.FileNames }
})
$btnTklFolder.Add_Click({
  $d = New-Object System.Windows.Forms.FolderBrowserDialog
  $d.Description = 'Select the folder holding the Kinase Library files'
  if ($d.ShowDialog() -eq 'OK') { Set-Tkl @($d.SelectedPath) }
})

# ---- Page 5: ready ----
$p5 = New-Page 'Ready to install'
$lblReady = Add-Label $p5 '' 0 330

# ---- Page 6: installing ----
$p6 = New-Page 'Installing'
$lblStep = Add-Label $p6 'Preparing...' 0 24
$bar = New-Object System.Windows.Forms.ProgressBar
$bar.Location = New-Object System.Drawing.Point(0, 30); $bar.Size = New-Object System.Drawing.Size(665, 18)
$bar.Style = 'Marquee'; $bar.MarqueeAnimationSpeed = 30
$p6.Controls.Add($bar)
$tbLog = New-Object System.Windows.Forms.TextBox
$tbLog.Multiline = $true; $tbLog.ReadOnly = $true; $tbLog.ScrollBars = 'Vertical'
$tbLog.Font = New-Object System.Drawing.Font('Consolas', 8.5)
$tbLog.Location = New-Object System.Drawing.Point(0, 58); $tbLog.Size = New-Object System.Drawing.Size(665, 320)
$p6.Controls.Add($tbLog)

# ---- Page 7: finish ----
$p7 = New-Page 'OSMSS is installed'
$tbSummary = New-Object System.Windows.Forms.TextBox
$tbSummary.Multiline = $true; $tbSummary.ReadOnly = $true; $tbSummary.ScrollBars = 'Vertical'
$tbSummary.Location = New-Object System.Drawing.Point(0, 0); $tbSummary.Size = New-Object System.Drawing.Size(665, 300)
$p7.Controls.Add($tbSummary)
$cbShortcut = New-Object System.Windows.Forms.CheckBox
$cbShortcut.Text = 'Create a desktop shortcut for OSMSS'; $cbShortcut.Checked = $true
$cbShortcut.Location = New-Object System.Drawing.Point(0, 310); $cbShortcut.Size = New-Object System.Drawing.Size(400, 24)
$cbStart = New-Object System.Windows.Forms.CheckBox
$cbStart.Text = 'Start OSMSS now'; $cbStart.Checked = $false
$cbStart.Location = New-Object System.Drawing.Point(0, 338); $cbStart.Size = New-Object System.Drawing.Size(400, 24)
$p7.Controls.Add($cbShortcut); $p7.Controls.Add($cbStart)

# ------------------------------------------------------------------------------------------
# Installation (runs on page 6)
# ------------------------------------------------------------------------------------------
$timer = New-Object System.Windows.Forms.Timer
$timer.Interval = 1000

function Show-Log {
  $lines = New-Object System.Collections.ArrayList
  foreach ($l in $state.InstallLog) { [void]$lines.Add($l) }
  foreach ($f in @($state.OutLog, $state.ErrLog)) {
    if ($f -and (Test-Path $f)) {
      try {
        $fs = [System.IO.File]::Open($f, 'Open', 'Read', 'ReadWrite')
        $sr = New-Object System.IO.StreamReader($fs)
        $txt = $sr.ReadToEnd(); $sr.Close(); $fs.Close()
        $tail = @($txt -split "`r?`n" | Where-Object { $_ -ne '' } | Select-Object -Last 12)
        foreach ($l in $tail) { [void]$lines.Add('  ' + $l) }
      } catch { }
    }
  }
  $tbLog.Text = ($lines -join "`r`n")
  $tbLog.SelectionStart = $tbLog.Text.Length; $tbLog.ScrollToCaret()
}

function Start-Install {
  $btnBack.Enabled = $false; $btnNext.Enabled = $false
  $t = $state.Target
  try {
    $lblStep.Text = 'Copying OSMSS files...'; Show-Log; [System.Windows.Forms.Application]::DoEvents()
    Copy-App $t
    if ($rbSaintNow.Checked -and $state.SaintDir) { $lblStep.Text = 'Setting up SAINTexpress...'; Install-Saint $t }
    if ($rbTklNow.Checked -and $state.Tkl.Count) { $lblStep.Text = 'Adding The Kinase Library files...'; Install-Tkl $t }
    $logDir = Join-Path $t 'logs'; New-Item -ItemType Directory -Force -Path $logDir | Out-Null
    $stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
    $state.OutLog = Join-Path $logDir "install_packages_$stamp.log"
    $state.ErrLog = Join-Path $logDir "install_packages_$stamp.err.log"
    $state.Report = Join-Path $logDir "install_report_$stamp.txt"
    $script = Join-Path $t 'installer\install_packages.R'
    if (-not (Test-Path $script)) { $script = Join-Path $SourceRoot 'installer\install_packages.R' }
    $lblStep.Text = 'Installing R packages (10 to 30 minutes; keep this window open)...'
    Add-Note "Running R $($state.R.Version): $script"
    $state.Started = Get-Date
    $state.Proc = Start-Process -FilePath $state.R.Rscript -ArgumentList @("`"$script`"", "`"$($state.Report)`"") `
      -WorkingDirectory $t -RedirectStandardOutput $state.OutLog -RedirectStandardError $state.ErrLog `
      -NoNewWindow -PassThru
    $null = $state.Proc.Handle   # keeps the exit code readable after the process ends
    $timer.Start()
  } catch {
    Add-Note ("ERROR: " + $_.Exception.Message); Show-Log
    $lblStep.Text = 'Installation stopped. See the messages below.'; $bar.Style = 'Blocks'
    $btnBack.Enabled = $true
  }
}

$timer.Add_Tick({
  Show-Log
  $mins = [int]((Get-Date) - $state.Started).TotalMinutes
  $lblStep.Text = "Installing R packages ($mins min so far; usually 10 to 30 minutes; keep this window open)..."
  if ($state.Proc.HasExited) {
    $timer.Stop(); $bar.Style = 'Blocks'; $bar.Value = 100
    Add-Note "R finished (exit code $($state.Proc.ExitCode))."
    $state.Done = $true
    Show-Finish
  }
})

function Show-Finish {
  $rep = Read-Report $state.Report
  $pre    = @(($rep['PREINSTALLED'] -split ',') | Where-Object { $_ })
  $new    = @(($rep['INSTALLED'] -split ',') | Where-Object { $_ })
  $failed = @(($rep['FAILED'] -split ',') | Where-Object { $_ })
  $s = New-Object System.Collections.ArrayList
  [void]$s.Add("Installed in: $($state.Target)")
  [void]$s.Add("R: $($state.R.Version)  ($($state.R.Rscript))")
  [void]$s.Add('')
  if (-not $rep.ContainsKey('STATUS')) {
    [void]$s.Add('The R package step did not finish. See the log:')
    [void]$s.Add("  $($state.OutLog)")
  } else {
    [void]$s.Add("R packages newly installed: $($new.Count)")
    if ($pre.Count) {
      [void]$s.Add("R packages you already had (left as they were, not updated): $($pre.Count)")
      [void]$s.Add('  ' + ($pre -join ', '))
      [void]$s.Add('  If OSMSS misbehaves, an old package version may be the reason; update packages in R with BiocManager::install().')
    }
    if ($failed.Count) {
      [void]$s.Add('')
      [void]$s.Add("PROBLEM: these packages could not be installed or loaded: $($failed -join ', ')")
      [void]$s.Add('  In R, try:  BiocManager::install(c("' + ($failed -join '", "') + '"))')
      [void]$s.Add("  Full log: $($state.OutLog)")
    } else {
      [void]$s.Add('All R packages OSMSS needs are installed.')
    }
  }
  [void]$s.Add('')
  [void]$s.Add('SAINTexpress (AP-MS): ' + $(if ($rbSaintNow.Checked -and $state.SaintDir) { 'set up' } else { 'not set up (can be added later; see README.md)' }))
  [void]$s.Add('The Kinase Library (kinase prediction): ' + $(if ($rbTklNow.Checked -and $state.Tkl.Count) { 'added' } else { 'not added (can be added later; see README.md)' }))
  [void]$s.Add('')
  [void]$s.Add('To start OSMSS, use the desktop shortcut or "Start OSMSS.cmd" in the install folder.')
  [void]$s.Add('OSMSS opens in your web browser; closing its black console window stops OSMSS.')
  $tbSummary.Text = ($s -join "`r`n")
  Show-Page 7
}

# ------------------------------------------------------------------------------------------
# Navigation
# ------------------------------------------------------------------------------------------
$script:page = 0
function Show-Page([int]$n) {
  for ($k = 0; $k -lt $pages.Count; $k++) { $pages[$k].Visible = ($k -eq $n) }
  $script:page = $n
  $header.Text = $pages[$n].Tag
  $btnBack.Enabled = ($n -gt 0 -and $n -lt 6)
  $btnNext.Enabled = $true
  $btnNext.Text = 'Next >'
  $btnCancel.Enabled = ($n -lt 7)
  switch ($n) {
    2 { if (-not $state.R) { $state.R = Find-R }; Update-RStatus }
    5 {
      $lines = @(
        "Install folder:  $($state.Target)",
        "R:  $($state.R.Version)",
        ('SAINTexpress:  ' + $(if ($rbSaintNow.Checked -and $state.SaintDir) { 'will be set up' } else { 'skipped' })),
        ('The Kinase Library:  ' + $(if ($rbTklNow.Checked -and $state.Tkl.Count) { 'will be added (' + ($state.Tkl.Keys -join ', ') + ')' } else { 'skipped' })),
        '',
        'Click Install to copy OSMSS and install the R packages it needs. Packages you already have are not updated.',
        'This needs an internet connection and usually takes 10 to 30 minutes.')
      $lblReady.Text = ($lines -join "`r`n`r`n")
      $btnNext.Text = 'Install'
    }
    6 { $btnNext.Enabled = $false; $btnBack.Enabled = $false; Start-Install }
    7 { $btnNext.Text = 'Finish'; $btnBack.Enabled = $false }
  }
}

$btnNext.Add_Click({
  switch ($script:page) {
    1 {
      $t = $tbFolder.Text.Trim()
      if (-not $t) { [void][System.Windows.Forms.MessageBox]::Show('Choose an install folder.', 'OSMSS Setup'); return }
      try { $t = [System.IO.Path]::GetFullPath($t) } catch { [void][System.Windows.Forms.MessageBox]::Show('That folder name is not valid.', 'OSMSS Setup'); return }
      if ((Test-Path $t) -and -not (Test-SamePath $t $SourceRoot)) {
        $items = @(Get-ChildItem $t -Force -ErrorAction SilentlyContinue)
        if ($items.Count -gt 0 -and -not (Test-Path (Join-Path $t 'app.R'))) {
          $ans = [System.Windows.Forms.MessageBox]::Show("The folder`r`n$t`r`nis not empty. Install OSMSS into it anyway?", 'OSMSS Setup', 'YesNo', 'Question')
          if ($ans -ne 'Yes') { return }
        }
      }
      try { New-Item -ItemType Directory -Force -Path $t | Out-Null } catch {
        [void][System.Windows.Forms.MessageBox]::Show("Cannot create or write to`r`n$t`r`nChoose another folder.", 'OSMSS Setup'); return
      }
      $state.Target = $t
    }
    3 { if ($rbSaintNow.Checked -and -not $state.SaintDir) { [void][System.Windows.Forms.MessageBox]::Show('Select the SAINTexpress download first, or choose "Skip this step".', 'OSMSS Setup'); return } }
    4 {
      $missing = @($TklRequired | Where-Object { -not $state.Tkl.ContainsKey($_) })
      if ($rbTklNow.Checked -and $missing.Count) { [void][System.Windows.Forms.MessageBox]::Show(('Still needed: ' + ($missing -join ', ') + "`r`nSelect them, or choose ""Skip this step""."), 'OSMSS Setup'); return }
    }
    7 {
      try {
        $launcher = Write-Launcher $state.Target
        if ($cbShortcut.Checked) { [void](New-DesktopShortcut $state.Target $launcher) }
        if ($cbStart.Checked) { Start-Process -FilePath $launcher -WorkingDirectory $state.Target }
      } catch { [void][System.Windows.Forms.MessageBox]::Show("Could not create the shortcut: $($_.Exception.Message)", 'OSMSS Setup') }
      $form.Close(); return
    }
  }
  Show-Page ($script:page + 1)
})
$btnBack.Add_Click({ if ($script:page -gt 0) { Show-Page ($script:page - 1) } })
$btnCancel.Add_Click({
  $msg = if ($state.Proc -and -not $state.Proc.HasExited) { 'R packages are still being installed. Stop and exit Setup?' } else { 'Exit OSMSS Setup?' }
  if ([System.Windows.Forms.MessageBox]::Show($msg, 'OSMSS Setup', 'YesNo', 'Question') -eq 'Yes') {
    if ($state.Proc -and -not $state.Proc.HasExited) { try { Stop-Process -Id $state.Proc.Id -Force } catch { } }
    $state.Proc = $null
    $form.Close()
  }
})
$form.Add_FormClosing({
  param($s, $e)
  if ($state.Proc -and -not $state.Proc.HasExited) {
    $ans = [System.Windows.Forms.MessageBox]::Show('R packages are still being installed. Stop and exit Setup?', 'OSMSS Setup', 'YesNo', 'Question')
    if ($ans -ne 'Yes') { $e.Cancel = $true; return }
    try { Stop-Process -Id $state.Proc.Id -Force } catch { }
  }
})
$form.Add_FormClosed({ $timer.Stop(); Remove-TempDirs })

# Files unpacked from a downloaded ZIP carry Windows' "downloaded from the internet" mark;
# clear it on the OSMSS files so R and the launcher run without extra prompts.
try {
  Get-ChildItem $SourceRoot -Force -ErrorAction SilentlyContinue |
    Where-Object { $ExcludeDirs -notcontains $_.Name } |
    ForEach-Object { if ($_.PSIsContainer) { Get-ChildItem $_.FullName -Recurse -File -ErrorAction SilentlyContinue } else { $_ } } |
    Unblock-File -ErrorAction SilentlyContinue
} catch { }

if ($UiTestReport) {
  # Test mode: build every page and the finish summary without showing the window.
  foreach ($n in 0..5) { Show-Page $n; "page ${n}: " + $pages[$n].Tag + " | Next enabled: " + $btnNext.Enabled }
  $state.Report = $UiTestReport; $state.OutLog = $UiTestReport
  Show-Finish
  "finish page: " + $header.Text
  $tbSummary.Text
  return
}

Show-Page 0
[void]$form.ShowDialog()
