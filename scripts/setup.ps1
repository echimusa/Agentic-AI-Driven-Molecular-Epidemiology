<#
.SYNOPSIS
  Install what the exercises of the chosen chapters need: a Python environment, R packages and, where available, the command-line tools.
.DESCRIPTION
  Ubuntu (or WSL2 Ubuntu): uses apt, a virtual environment in .venv and a user R library in R/library. Needs sudo for apt.
  Windows: uses winget for Python 3.12, R, Pandoc and Git, then the same virtual environment and R library, for the chapters that run natively.
  Chapters that need Linux tools (Chapter 5: IQ-TREE 2 and FastTree; Chapters 6 and 7: mappy; with -Rebuild or -Heavy also MAFFT, bcftools, PLINK, BEAST 2 and
  AMRFinderPlus) are prepared inside WSL2 Ubuntu: the script finds the Ubuntu distribution, copies the repository into its home folder and runs scripts/setup.sh there.
  If WSL2 is not installed it offers to install it (administrator approval and a restart are needed; then run the script again).
.PARAMETER Chapters  Chapters to prepare (default 1 to 13).
.PARAMETER Heavy     Also prepare the heavy analyses: check BEAST 2 (it brings its own Java), the MODEL_SELECTION package, AMRFinderPlus and its database, and download and install what is missing (scripts/setup-heavy.sh; on Windows inside WSL2).
.PARAMETER Rebuild   Also install the tools needed to rebuild alignments, trees and genotype extracts from public sources.
.PARAMETER SkipSystem  Do not touch apt or winget (use this when Python, R and the tools are already installed).
.PARAMETER SkipWsl   Windows only: do not use WSL2. The chapters that need Linux tools are then reported as not ready.
.PARAMETER InstallWsl  Windows only: install WSL2 with Ubuntu without asking, if no Ubuntu distribution is found.
.PARAMETER DryRun    Print every step and run none.
.EXAMPLE
  pwsh scripts/setup.ps1 -Chapters 12,13
  pwsh scripts/setup.ps1 -Chapters 5,6,7 -InstallWsl
  pwsh scripts/setup.ps1 -Rebuild
#>
[CmdletBinding()]
param([string[]]$Chapters = @('1-13'), [switch]$Heavy, [switch]$Rebuild, [switch]$SkipSystem, [switch]$SkipWsl, [switch]$InstallWsl, [switch]$DryRun)
. (Join-Path $PSScriptRoot '_common.ps1')
$Chapters = ConvertTo-ChapterList $Chapters
$isWin = Test-IsWindows

# On Windows, chapters that need Linux tools go to WSL2; the others are prepared natively
$wslChapters = @()
if ($isWin -and -not $SkipWsl) { $wslChapters = @(Get-WslChapters $Chapters -Rebuild:$Rebuild -Heavy:$Heavy) }
$nativeChapters = @($Chapters | Where-Object { $_ -notin $wslChapters })
Write-Host "Chapters: $($Chapters -join ', ')"
if ($wslChapters.Count) {
    Write-Host "  prepared on Windows: $(if ($nativeChapters.Count) { $nativeChapters -join ', ' } else { '(none)' })"
    Write-Host "  prepared in WSL2 Ubuntu (Linux tools needed): $($wslChapters -join ', ')"
}

if ($nativeChapters.Count) {
    $need = Get-NeededItems $nativeChapters
    Write-Host "Python packages: $($need.Python -join ', ')"
    Write-Host "R packages: $(if ($need.R.Count) { $need.R -join ', ' } else { '(none: base R only or no R)' })"
    Write-Host "Tools needed for the quick check: $(if ($need.Quick.Count) { $need.Quick -join ', ' } else { '(none)' })"

    if (-not $SkipSystem) {
        if ($isWin) {
            if (-not (Get-Command winget -ErrorAction SilentlyContinue)) { throw 'winget was not found. Install "App Installer" from the Microsoft Store, or install Python 3.12, R and Pandoc by hand and rerun with -SkipSystem.' }
            $spec = Get-Spec
            $ids = @($spec.winget.python, $spec.winget.pandoc, $spec.winget.git)
            if ($need.R.Count -or $nativeChapters.Count -gt 0) { $ids += $spec.winget.r }
            if ($Heavy -and -not $wslChapters.Count) { $ids += $spec.winget.java }
            foreach ($id in $ids) { Invoke-Step "winget install $id" { winget install --id $id -e --accept-package-agreements --accept-source-agreements } -DryRun:$DryRun }
            if ($SkipWsl) {
                if ($need.Quick.Count) { Write-Warning "Chapter(s) needing $($need.Quick -join ', ') have no native Windows route here. Run without -SkipWsl to prepare them in WSL2 (docs/INSTALL-Windows.md)." }
                if ($Rebuild -or $need.Rebuild.Count) { Write-Warning 'Rebuilding from public sources needs MAFFT, FastTree, bcftools, tabix or PLINK: use WSL2. The quick checks do not need them.' }
                if ('mappy' -in $need.Python) { Write-Warning 'mappy (Chapters 6 and 7) has no Windows wheel and needs the Visual C++ Build Tools to compile. If pip fails for it, use WSL2 for Chapters 6 and 7.' }
            }
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

    # BEAST 2, Java, AMRFinderPlus and its database: check each one, install what is missing (Ubuntu; on Windows these chapters are prepared in WSL2)
    if ($Heavy -and -not $isWin) {
        Write-Host "==> check BEAST 2, Java and AMRFinderPlus, and install what is missing" -ForegroundColor Cyan
        $hargs = @('--chapters', ($nativeChapters -join ',')); if ($DryRun) { $hargs += '--dry-run' }
        bash (Join-Path $PSScriptRoot 'setup-heavy.sh') @hargs
        if ($LASTEXITCODE -ne 0) { Write-Warning 'BEAST 2 or AMRFinderPlus is not fully installed: read the messages above, fix the cause and run scripts/setup-heavy.sh again.' }
    }
}

# Chapters that need Linux tools: prepare them in WSL2 Ubuntu
$wslDistro = $null; $wslDone = $false; $wslInstalled = $false
if ($wslChapters.Count) {
    Write-Host "`n==> Chapter(s) $($wslChapters -join ', ') need Linux tools: preparing them in WSL2 Ubuntu" -ForegroundColor Cyan
    $wslDistro = Get-WslDistro
    if (-not $wslDistro) {
        $install = [bool]$InstallWsl
        if (-not $install -and -not $DryRun) {
            try { $install = ((Read-Host 'WSL2 with Ubuntu was not found. Install WSL2 and Ubuntu now? It needs administrator approval and a restart [y/N]') -match '^(y|yes|j|o|oui)$') } catch { $install = $false }
        }
        if ($install -or $DryRun) {
            Invoke-Step 'wsl --install -d Ubuntu --no-launch   (asks for administrator approval)' {
                # 3010 means "done, restart needed". If this WSL version does not know --no-launch, the second attempt installs and launches Ubuntu.
                $ok = $false
                foreach ($extra in @(@('--no-launch'), @())) {
                    $p = Start-Process -FilePath 'wsl.exe' -ArgumentList (@('--install', '-d', 'Ubuntu') + $extra) -Verb RunAs -Wait -PassThru
                    if ($p.ExitCode -in 0, 3010) { $ok = $true; break }
                }
                $script:wslInstalled = $ok
                if (-not $ok) { Write-Warning "wsl --install did not finish (code $($p.ExitCode)). Run 'wsl --install -d Ubuntu' by hand in an administrator PowerShell; 'wsl --list --online' shows the names this PC accepts (docs/INSTALL-Windows.md, Route A)." }
            } -DryRun:$DryRun
            if ($wslInstalled -or $DryRun) { Write-Host 'Next: restart Windows if it asks, open "Ubuntu" once from the Start menu to create your Linux user, then run this script again.' -ForegroundColor Yellow }
        } else {
            Write-Warning "Chapter(s) $($wslChapters -join ', ') were not prepared: WSL2 with Ubuntu is not installed. Run 'wsl --install -d Ubuntu' in an administrator PowerShell, restart, open Ubuntu once, then run this script again (or add -InstallWsl)."
        }
    } else {
        Invoke-Step "set up chapter(s) $($wslChapters -join ', ') in $wslDistro" {
            $rc = Invoke-WslSetup -Distro $wslDistro -Chapters $wslChapters -Rebuild:$Rebuild -Heavy:$Heavy
            if ($rc -ne 0) { Write-Warning "The setup in WSL2 ended with code ${rc}: see the messages above and docs/TROUBLESHOOTING.md." } else { $script:wslDone = $true }
        } -DryRun:$DryRun
    }
}
if (-not $DryRun) {
    if ($nativeChapters.Count) {
        Write-Host "`n==> check the installation" -ForegroundColor Cyan
        & (Join-Path $PSScriptRoot 'check-install.ps1') -Chapters $nativeChapters
    }
    if ($wslDone) {
        Write-Host "`nChapter(s) $($wslChapters -join ', ') run inside WSL2: open Ubuntu, cd ~/$(Split-Path $RepoRoot -Leaf) and use  bash scripts/run-chapter.sh <number>" -ForegroundColor Green
    }
}
