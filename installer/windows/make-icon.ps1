# Builds installer\windows\OSMSS.ico (the desktop-shortcut icon) in the app's colours:
# a navy rounded square with a mass-spectrum peak plot. Maintainer tool; the .ico is committed.
# Usage:  powershell -NoProfile -ExecutionPolicy Bypass -File make-icon.ps1 [-PreviewPng path]
param([string]$PreviewPng)
Add-Type -AssemblyName System.Drawing

$navy = [System.Drawing.Color]::FromArgb(255, 20, 61, 76)     # #143d4c
$teal = [System.Drawing.Color]::FromArgb(255, 24, 188, 156)   # #18bc9c
$pale = [System.Drawing.Color]::FromArgb(255, 220, 240, 236)

function New-IconBitmap([int]$s) {
  $bmp = New-Object System.Drawing.Bitmap($s, $s, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
  $g = [System.Drawing.Graphics]::FromImage($bmp)
  $g.SmoothingMode = 'AntiAlias'; $g.PixelOffsetMode = 'HighQuality'
  $g.Clear([System.Drawing.Color]::Transparent)
  # rounded square
  $r = [Math]::Max(2, [int]($s * 0.18)); $pad = [Math]::Max(0, [int]($s * 0.03)); $w = $s - 2 * $pad
  $path = New-Object System.Drawing.Drawing2D.GraphicsPath
  $path.AddArc($pad, $pad, 2 * $r, 2 * $r, 180, 90)
  $path.AddArc($pad + $w - 2 * $r, $pad, 2 * $r, 2 * $r, 270, 90)
  $path.AddArc($pad + $w - 2 * $r, $pad + $w - 2 * $r, 2 * $r, 2 * $r, 0, 90)
  $path.AddArc($pad, $pad + $w - 2 * $r, 2 * $r, 2 * $r, 90, 90)
  $path.CloseFigure()
  $g.FillPath((New-Object System.Drawing.SolidBrush($navy)), $path)
  # spectrum: baseline and peaks (relative heights), the tallest in teal
  $left = $s * 0.20; $right = $s * 0.80; $base = $s * 0.76; $top = $s * 0.20
  $heights = if ($s -le 16) { @(0.45, 1.0, 0.6) } elseif ($s -le 32) { @(0.35, 0.7, 1.0, 0.5) } else { @(0.3, 0.55, 1.0, 0.4, 0.7, 0.25) }
  $n = $heights.Count; $step = ($right - $left) / ($n - 1)
  $pw = [Math]::Max(1.6, $s * $(if ($s -le 16) { 0.13 } elseif ($s -le 32) { 0.09 } else { 0.065 }))
  $tallest = [Array]::IndexOf($heights, ($heights | Measure-Object -Maximum).Maximum)
  for ($k = 0; $k -lt $n; $k++) {
    $x = $left + $k * $step; $y = $base - $heights[$k] * ($base - $top)
    $col = if ($k -eq $tallest) { $teal } else { $pale }
    $pen = New-Object System.Drawing.Pen($col, $pw); $pen.StartCap = 'Flat'; $pen.EndCap = 'Flat'
    $g.DrawLine($pen, [single]$x, [single]$base, [single]$x, [single]$y)
    $pen.Dispose()
  }
  $bpen = New-Object System.Drawing.Pen($pale, [Math]::Max(1.2, $s * 0.035))
  $g.DrawLine($bpen, [single]($s * 0.14), [single]($base + $pw * 0.1), [single]($s * 0.86), [single]($base + $pw * 0.1))
  $g.Dispose()
  return $bmp
}

$sizes = @(256, 64, 48, 32, 24, 16)
$pngs = foreach ($sz in $sizes) {
  $b = New-IconBitmap $sz; $ms = New-Object System.IO.MemoryStream
  $b.Save($ms, [System.Drawing.Imaging.ImageFormat]::Png); $b.Dispose()
  ,$ms.ToArray()
}
# ICO container with PNG-compressed images (supported by Windows Vista and later)
$out = Join-Path $PSScriptRoot 'OSMSS.ico'
$fs = [System.IO.File]::Create($out); $bw = New-Object System.IO.BinaryWriter($fs)
$bw.Write([uint16]0); $bw.Write([uint16]1); $bw.Write([uint16]$sizes.Count)
$offset = 6 + 16 * $sizes.Count
for ($k = 0; $k -lt $sizes.Count; $k++) {
  $d = if ($sizes[$k] -ge 256) { 0 } else { $sizes[$k] }
  $bw.Write([byte]$d); $bw.Write([byte]$d); $bw.Write([byte]0); $bw.Write([byte]0)
  $bw.Write([uint16]1); $bw.Write([uint16]32); $bw.Write([uint32]$pngs[$k].Length); $bw.Write([uint32]$offset)
  $offset += $pngs[$k].Length
}
foreach ($p in $pngs) { $bw.Write($p) }
$bw.Close(); $fs.Close()
"wrote $out ($((Get-Item $out).Length) bytes)"

# macOS (.icns) and Linux (.png) icons, in installer\icons
$iconDir = Join-Path (Split-Path $PSScriptRoot) 'icons'
New-Item -ItemType Directory -Force -Path $iconDir | Out-Null
function Get-PngBytes([int]$sz) {
  $b = New-IconBitmap $sz; $ms = New-Object System.IO.MemoryStream
  $b.Save($ms, [System.Drawing.Imaging.ImageFormat]::Png); $b.Dispose(); return ,$ms.ToArray()
}
$png = Join-Path $iconDir 'OSMSS.png'
[System.IO.File]::WriteAllBytes($png, (Get-PngBytes 256))
"wrote $png"
# ICNS: 'icns' + big-endian total length, then entries of type (4 chars) + big-endian length + PNG data
function Write-BE([System.IO.BinaryWriter]$w, [uint32]$v) { $bytes = [BitConverter]::GetBytes($v); [Array]::Reverse($bytes); $w.Write($bytes) }
$entries = @(@('icp4', 16), @('icp5', 32), @('icp6', 64), @('ic07', 128), @('ic08', 256), @('ic09', 512), @('ic10', 1024))
$data = foreach ($e in $entries) { ,@($e[0], (Get-PngBytes $e[1])) }
$total = 8; foreach ($d in $data) { $total += 8 + $d[1].Length }
$icns = Join-Path $iconDir 'OSMSS.icns'
$fs = [System.IO.File]::Create($icns); $bw = New-Object System.IO.BinaryWriter($fs)
$bw.Write([System.Text.Encoding]::ASCII.GetBytes('icns')); Write-BE $bw $total
foreach ($d in $data) { $bw.Write([System.Text.Encoding]::ASCII.GetBytes($d[0])); Write-BE $bw (8 + $d[1].Length); $bw.Write($d[1]) }
$bw.Close(); $fs.Close()
"wrote $icns ($((Get-Item $icns).Length) bytes)"

if ($PreviewPng) {
  # one strip with every size, for checking the design
  $strip = New-Object System.Drawing.Bitmap(460, 270)
  $g = [System.Drawing.Graphics]::FromImage($strip); $g.Clear([System.Drawing.Color]::White)
  $b = New-IconBitmap 256; $g.DrawImage($b, 5, 5, 256, 256); $b.Dispose()
  $x = 275; $y = 5
  foreach ($sz in ($sizes | Where-Object { $_ -lt 256 })) {
    $b = New-IconBitmap $sz; $g.DrawImage($b, $x, $y, $sz, $sz); $b.Dispose()
    $x += $sz + 10; if ($x -gt 400) { $x = 275; $y += 70 }
  }
  $g.Dispose(); $strip.Save($PreviewPng, [System.Drawing.Imaging.ImageFormat]::Png); $strip.Dispose()
  "preview $PreviewPng"
}
