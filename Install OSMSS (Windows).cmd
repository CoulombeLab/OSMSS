@echo off
rem Double-click this file to install OSMSS on Windows.
rem It opens the OSMSS Setup wizard (installer\windows\Install-OSMSS.ps1).
start "OSMSS Setup" /min powershell.exe -NoProfile -ExecutionPolicy Bypass -STA -WindowStyle Hidden -File "%~dp0installer\windows\Install-OSMSS.ps1"
