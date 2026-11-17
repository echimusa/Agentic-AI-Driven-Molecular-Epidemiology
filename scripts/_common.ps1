# Shared helpers, dot-sourced by the other scripts. Works in Windows PowerShell 5.1 and in PowerShell 7 on Windows or Ubuntu.
$ErrorActionPreference = 'Stop'
$script:RepoRoot = Split-Path -Parent $PSScriptRoot
function Test-IsWindows { if ($PSVersionTable.PSEdition -eq 'Desktop') { return $true }; return [bool]$IsWindows }
function Get-Spec { Get-Content (Join-Path $script:RepoRoot 'tools/chapter_requirements.json') -Raw | ConvertFrom-Json }
function Get-VenvPython {
    if (Test-IsWindows) { Join-Path $script:RepoRoot '.venv\Scripts\python.exe' } else { Join-Path $script:RepoRoot '.venv/bin/python' }
}
function Get-Python {
    $v = Get-VenvPython
    if (Test-Path $v) { return $v }
    foreach ($c in @('python3', 'python')) {
        $cmd = Get-Command $c -ErrorAction SilentlyContinue
        if ($cmd -and $cmd.Source -notmatch 'WindowsApps') { return $cmd.Source }
    }
    return $null
}
function Get-Rscript {
    $cmd = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    if (Test-IsWindows) {
        $f = Get-ChildItem 'C:\Program Files\R\R-*\bin\Rscript.exe' -ErrorAction SilentlyContinue | Sort-Object FullName -Descending | Select-Object -First 1
        if ($f) { return $f.FullName }
    }
    return $null
}
function Set-RunEnvironment {
    # R packages go to a user library inside the repository folder, so nothing needs administrator rights
    $lib = Join-Path $script:RepoRoot 'R/library'
    New-Item -ItemType Directory -Force -Path $lib | Out-Null
    $env:R_LIBS_USER = $lib
    $env:PYTHONUTF8 = '1'
    $env:MPLBACKEND = 'Agg'
    $sep = [IO.Path]::PathSeparator
    $extra = @((Join-Path $script:RepoRoot 'tools/linux/bin'), (Join-Path $script:RepoRoot 'tools/win'))
    foreach ($e in $extra) { if ((Test-Path $e) -and ($env:PATH -notlike "*$e*")) { $env:PATH = "$e$sep$env:PATH" } }
    $r = Get-Rscript
    if ($r) { $d = Split-Path -Parent $r; if ($env:PATH -notlike "*$d*") { $env:PATH = "$d$sep$env:PATH" } }
}
function Get-NeededItems([int[]]$Chapters) {
    $spec = Get-Spec
    $py = [System.Collections.Generic.HashSet[string]]::new(); $r = [System.Collections.Generic.HashSet[string]]::new()
    $quick = [System.Collections.Generic.HashSet[string]]::new(); $rebuild = [System.Collections.Generic.HashSet[string]]::new(); $heavy = [System.Collections.Generic.HashSet[string]]::new()
    foreach ($p in $spec.base.python) { [void]$py.Add($p) }
    foreach ($c in $Chapters) {
        $e = $spec.chapters.PSObject.Properties[[string]$c].Value
        foreach ($p in $e.python) { [void]$py.Add($p) }
        foreach ($p in $e.r) { [void]$r.Add($p) }
        foreach ($p in $e.quick_tools) { [void]$quick.Add($p) }
        foreach ($p in $e.rebuild_tools) { [void]$rebuild.Add($p) }
        if ($e.PSObject.Properties['heavy_tools']) { foreach ($p in $e.heavy_tools) { [void]$heavy.Add($p) } }
    }
    [pscustomobject]@{ Python = @($py); R = @($r); Quick = @($quick); Rebuild = @($rebuild); Heavy = @($heavy) }
}
function Invoke-Step([string]$Text, [scriptblock]$Action, [switch]$DryRun) {
    Write-Host "==> $Text" -ForegroundColor Cyan
    if ($DryRun) { Write-Host "    (dry run: not executed)" -ForegroundColor DarkGray; return }
    & $Action
}
function ConvertTo-ChapterList($Value) {
    # accepts 12,13  or  "12,13"  or  1-5  or  1..5  (pwsh -File passes a comma list as one string)
    $out = @()
    foreach ($t in @($Value)) {
        foreach ($p in ([string]$t -split '[,\s]+')) {
            if (-not $p) { continue }
            if ($p -match '^(\d+)(\.\.|-)(\d+)$') { $out += [int]$Matches[1]..[int]$Matches[3] } else { $out += [int]$p }
        }
    }
    @($out | Sort-Object -Unique)
}
