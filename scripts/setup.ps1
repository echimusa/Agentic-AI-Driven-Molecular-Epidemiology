<#
.SYNOPSIS
  Install what the exercises of the chosen chapters need: a Python environment, R packages and, where available, the command-line tools.
.DESCRIPTION
  Ubuntu (or WSL2 Ubuntu): uses apt, a virtual environment in .venv and a user R library in R/library. Needs sudo for apt.
  Windows: uses winget for Python 3.12, R, Pandoc and Git, then the same virtual environment and R library.
  Command-line tools that have no Windows build (IQ-TREE 2 here is the one that matters, for Chapter 5; also MAFFT, FastTree, bcftools, PLINK) are not installed
  on native Windows; use WSL2 for those chapters (docs/INSTALL-Windows.md).
.PARAMETER Chapters  Chapters to prepare (default 1 to 13).
.PARAMETER Heavy     Also prepare the heavy analyses (Java for BEAST 2, BLAST and HMMER for AMRFinderPlus); then run scripts/setup-heavy.ps1.
.PARAMETER Rebuild   Also install the tools needed to rebuild alignments, trees and genotype extracts from public sources.
.PARAMETER SkipSystem  Do not touch apt or winget (use this when Python, R and the tools are already installed).
.PARAMETER DryRun    Print every step and run none.
.EXAMPLE
  pwsh scripts/setup.ps1 -Chapters 12,13
  pwsh scripts/setup.ps1 -Rebuild
#>
[CmdletBinding()]
param([string[]]$Chapters = @('1-13'), [switch]$Heavy, [switch]$Rebuild, [switch]$SkipSystem, [switch]$DryRun)
. (Join-Path $PSScriptRoot '_common.ps1')
$Chapters = ConvertTo-ChapterList $Chapters
$need = Get-NeededItems $Chapters
$isWin = Test-IsWindows
Write-Host "Chapters: $($Chapters -join ', ')"
Write-Host "Python packages: $($need.Python -join ', ')"
Write-Host "R packages: $(if ($need.R.Count) { $need.R -join ', ' } else { '(none: base R only or no R)' })"
Write-Host "Tools needed for the quick check: $(if ($need.Quick.Count) { $need.Quick -join ', ' } else { '(none)' })"

if (-not $SkipSystem) {
    if ($isWin) {
        if (-not (Get-Command winget -ErrorAction SilentlyContinue)) { throw 'winget was not found. Install "App Installer" from the Microsoft Store, or install Python 3.12, R and Pandoc by hand and rerun with -SkipSystem.' }
        $spec = Get-Spec
        $ids = @($spec.winget.python, $spec.winget.pandoc, $spec.winget.git)
        if ($need.R.Count -or $Chapters.Count -gt 0) { $ids += $spec.winget.r }
        if ($Heavy) { $ids += $spec.winget.java }
        foreach ($id in $ids) { Invoke-Step "winget install $id" { winget install --id $id -e --accept-package-agreements --accept-source-agreements } -DryRun:$DryRun }
        if ($need.Quick.Count) { Write-Warning "Chapter(s) needing $($need.Quick -join ', ') have no native Windows route here. Use WSL2 for them (docs/INSTALL-Windows.md)." }
        if ($Rebuild -or $need.Rebuild.Count) { Write-Warning 'Rebuilding from public sources needs MAFFT, FastTree, bcftools, tabix or PLINK: use WSL2. The quick checks do not need them.' }
        if ('mappy' -in $need.Python) { Write-Warning 'mappy (Chapters 6 and 7) has no Windows wheel and needs the Visual C++ Build Tools to compile. If pip fails for it, use WSL2 for Chapters 6 and 7.' }
    } else {
        $sudo = if ((id -u) -eq '0') { '' } else { 'sudo ' }
        $apt = [System.Collections.Generic.List[string]]::new()
        'python3', 'python3-venv', 'python3-pip', 'python3-dev', 'build-essential', 'zlib1g-dev', 'r-base', 'r-base-dev', 'pandoc', 'curl', 'unzip', 'git', 'libcurl4-openssl-dev', 'libssl-dev', 'libxml2-dev', 'libfontconfig1-dev', 'libfreetype6-dev', 'libharfbuzz-dev', 'libfribidi-dev', 'libpng-dev', 'libjpeg-dev' | ForEach-Object { $apt.Add($_) }
        $spec = Get-Spec
        $tools = @($need.Quick)
        if ($Rebuild) { $tools += @($need.Rebuild) }
        if ($Heavy) { $tools += @($need.Heavy) }
        foreach ($t in ($tools | Select-Object -Unique)) { if ($spec.apt.PSObject.Properties[$t]) { foreach ($p in ($spec.apt.$t -split ' ')) { if ($p -and -not $apt.Contains($p)) { $apt.Add($p) } } } }
        Invoke-Step "apt-get update" { bash -c "${sudo}apt-get update" } -DryRun:$DryRun
        Invoke-Step "apt-get install $($apt -join ' ')" { bash -c "${sudo}DEBIAN_FRONTEND=noninteractive apt-get install -y $($apt -join ' ')" } -DryRun:$DryRun
    }
}

# Python environment
$pins = @{}
Get-Content (Join-Path $RepoRoot 'requirements.txt') | Where-Object { $_ -match '==' } | ForEach-Object { $a = $_.Trim() -split '=='; $pins[$a[0].ToLower()] = $_.Trim() }
$pipList = @($need.Python | ForEach-Object { if ($pins.ContainsKey($_.ToLower())) { $pins[$_.ToLower()] } else { $_ } })
Invoke-Step "create the virtual environment in .venv" {
    $py = Get-Python
    if (-not $py) { throw 'Python was not found. Open a new PowerShell window after installing it, or install Python 3.12 and rerun.' }
    if (-not (Test-Path (Get-VenvPython))) { & $py -m venv (Join-Path $RepoRoot '.venv') }
    & (Get-VenvPython) -m pip install --upgrade pip
} -DryRun:$DryRun
Invoke-Step "pip install $($pipList -join ' ')" {
    $failed = @()
    foreach ($p in $pipList) {
        & (Get-VenvPython) -m pip install $p
        if ($LASTEXITCODE -ne 0) { $failed += $p }
    }
    if ($failed.Count) { Write-Warning "pip could not install: $($failed -join ', '). See docs/TROUBLESHOOTING.md." }
} -DryRun:$DryRun

# R packages
if ($need.R.Count) {
    Invoke-Step "install R packages: $($need.R -join ', ')" {
        Set-RunEnvironment
        $rs = Get-Rscript
        if (-not $rs) { throw 'Rscript was not found. Install R (apt install r-base, or winget install RProject.R) and open a new PowerShell window.' }
        & $rs (Join-Path $RepoRoot 'r-packages.R') @($need.R)
        if ($LASTEXITCODE -ne 0) { Write-Warning 'Some R packages did not install; see the lines marked MISSING above and docs/TROUBLESHOOTING.md.' }
    } -DryRun:$DryRun
}
if ($Heavy) { Write-Host 'Next, for BEAST 2 and AMRFinderPlus: pwsh scripts/setup-heavy.ps1 (Ubuntu or WSL2 only).' -ForegroundColor Yellow }
if (-not $DryRun) {
    Write-Host "`n==> check the installation" -ForegroundColor Cyan
    & (Join-Path $PSScriptRoot 'check-install.ps1') -Chapters $Chapters
}
