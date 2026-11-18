<#
.SYNOPSIS  Check BEAST 2 (it brings its own Java), its MODEL_SELECTION package, AMRFinderPlus and its database, and download and install whatever is missing. Ubuntu or WSL2 only.
.DESCRIPTION
  Runs scripts/setup-heavy.sh, which does the work: it finds out what is already installed (nothing is downloaded twice), downloads only the missing parts, and tests that each one starts.
  On Windows these Linux builds are prepared inside WSL2 by:  .\scripts\setup.ps1 -Chapters 6,9 -Heavy
  RGI (Chapter 9, Hands-on 6) is not on PyPI: see docs/HEAVY-ANALYSES.md.
.PARAMETER Chapters  Prepare only what these chapters need (6: BEAST 2 and Java; 9: AMRFinderPlus). Default: both.
.PARAMETER Force     Download and install again, even if the tools are found.
.EXAMPLE   pwsh scripts/setup-heavy.ps1 -Chapters 9
#>
[CmdletBinding()]
param([string[]]$Chapters, [switch]$SkipBeast, [switch]$SkipAmrfinder, [switch]$Force, [switch]$DryRun)
. (Join-Path $PSScriptRoot '_common.ps1')
if (Test-IsWindows) { throw 'BEAST 2 and AMRFinderPlus are installed here as Linux builds. Run  .\scripts\setup.ps1 -Chapters 6,9 -Heavy  to prepare them inside WSL2 Ubuntu (docs/INSTALL-Windows.md).' }
$a = @()
if ($Chapters) { $a += '--chapters'; $a += (ConvertTo-ChapterList $Chapters) -join ',' }
if ($SkipBeast) { $a += '--skip-beast' }
if ($SkipAmrfinder) { $a += '--skip-amrfinder' }
if ($Force) { $a += '--force' }
if ($DryRun) { $a += '--dry-run' }
bash (Join-Path $PSScriptRoot 'setup-heavy.sh') @a
exit $LASTEXITCODE
