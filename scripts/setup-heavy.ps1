<#
.SYNOPSIS  Download BEAST 2 (with the model-selection package) and AMRFinderPlus for the heavy analyses of Chapters 6 and 9. Ubuntu or WSL2 only.
.NOTES     BEAST 2 needs Java (setup.ps1 -Heavy installs a headless JRE). The AMRFinderPlus database is large and changes with every release; the chapter used 4.2/2026-10-29.1.
           RGI (Chapter 9, Hands-on 6) is not on PyPI: see docs/HEAVY-ANALYSES.md.
#>
[CmdletBinding()]
param([switch]$DryRun, [switch]$SkipBeast, [switch]$SkipAmrfinder)
. (Join-Path $PSScriptRoot '_common.ps1')
if (Test-IsWindows) { throw 'These downloads are Linux builds. Run this script inside WSL2 Ubuntu (docs/INSTALL-Windows.md).' }
$bin = Join-Path $RepoRoot 'tools/linux'; New-Item -ItemType Directory -Force -Path $bin | Out-Null
if (-not $SkipBeast) {
    Invoke-Step 'download and unpack BEAST 2.7.7' {
        $tgz = Join-Path $bin 'beast.tgz'
        curl -L -o $tgz https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.x86.tgz
        tar -xzf $tgz -C $bin; Remove-Item $tgz
        New-Item -ItemType Directory -Force -Path (Join-Path $bin 'bin') | Out-Null
        foreach ($t in 'beast', 'beauti', 'treeannotator', 'packagemanager', 'loganalyser') { $src = Join-Path $bin "beast/bin/$t"; if (Test-Path $src) { Copy-Item $src (Join-Path $bin "bin/$t") -Force; (Get-Item (Join-Path $bin "bin/$t")).UnixFileMode = 'UserRead,UserWrite,UserExecute,GroupRead,GroupExecute,OtherRead,OtherExecute' } }
    } -DryRun:$DryRun
    Invoke-Step 'install the MODEL_SELECTION package 1.6.2 into ~/.beast/2.7' {
        $d = Join-Path $HOME '.beast/2.7/MODEL_SELECTION'; New-Item -ItemType Directory -Force -Path $d | Out-Null
        $z = Join-Path $bin 'ms.zip'; curl -L -o $z https://github.com/BEAST2-Dev/model-selection/releases/download/v1.6.0/MODEL_SELECTION.addon.v1.6.2.zip
        unzip -o $z -d $d | Out-Null; Remove-Item $z
    } -DryRun:$DryRun
}
if (-not $SkipAmrfinder) {
    Invoke-Step 'download AMRFinderPlus 4.2.7 binaries' {
        $a = Join-Path $RepoRoot 'amrf'; New-Item -ItemType Directory -Force -Path $a | Out-Null
        $t = Join-Path $a 'amrf.tgz'; curl -L -o $t https://github.com/ncbi/amr/releases/download/amrfinder_v4.2.7/amrfinder_binaries_v4.2.7.tar.gz
        tar -xzf $t -C $a; Remove-Item $t
    } -DryRun:$DryRun
    Write-Host 'Then fetch the database files (docs/HEAVY-ANALYSES.md lists the 43 files of release 2026-10-29.1) into amrf/db and run: amrf/amrfinder_index amrf/db' -ForegroundColor Yellow
}
Write-Host "Add $bin/bin to PATH for this session with: `$env:PATH = '$bin/bin:' + `$env:PATH"
