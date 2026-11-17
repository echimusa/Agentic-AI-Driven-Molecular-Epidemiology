<#
.SYNOPSIS  Report what is installed and which chapters are ready to run from the shipped data.
.EXAMPLE   pwsh scripts/check-install.ps1 -Chapters 12,13
#>
[CmdletBinding()]
param([string[]]$Chapters = @('1-13'))
. (Join-Path $PSScriptRoot '_common.ps1')
$Chapters = ConvertTo-ChapterList $Chapters
Set-RunEnvironment
$spec = Get-Spec; $need = Get-NeededItems $Chapters
$pins = @{}
Get-Content (Join-Path $RepoRoot 'requirements.txt') | Where-Object { $_ -match '==' } | ForEach-Object { $a = $_.Trim() -split '=='; $pins[$a[0].ToLower()] = $a[1] }
$missing = [System.Collections.Generic.HashSet[string]]::new(); $rows = @()
$py = Get-Python
if (-not $py) { $rows += [pscustomobject]@{ Item = 'python 3.12'; Status = 'MISSING'; Found = ''; Tested = '3.12.3' }; [void]$missing.Add('python') }
else {
    $v = (& $py -c 'import platform; print(platform.python_version())').Trim()
    $rows += [pscustomobject]@{ Item = 'python'; Status = $(if ($v -like '3.1[2-9]*') { 'OK' } else { 'OLD' }); Found = $v; Tested = '3.12.3' }
    foreach ($p in $need.Python) {
        $f = (& $py -c "import importlib.metadata as m,sys
try: print(m.version('$p'))
except Exception: print('')" 2>$null | Out-String).Trim()
        $t = $pins[$p.ToLower()]
        $st = if (-not $f) { 'MISSING' } elseif ($t -and $f -ne $t) { 'DIFFERENT' } else { 'OK' }
        if ($st -eq 'MISSING') { [void]$missing.Add($p) }
        $rows += [pscustomobject]@{ Item = "py: $p"; Status = $st; Found = $f; Tested = $t }
    }
}
$rs = Get-Rscript
if ($need.R.Count) {
    if (-not $rs) { $rows += [pscustomobject]@{ Item = 'R'; Status = 'MISSING'; Found = ''; Tested = '4.3.3' }; foreach ($p in $need.R) { [void]$missing.Add($p) } }
    else {
        $rv = (& $rs -e 'cat(as.character(getRversion()))' 2>$null | Out-String).Trim()
        $rows += [pscustomobject]@{ Item = 'R'; Status = 'OK'; Found = $rv; Tested = '4.3.3' }
        foreach ($p in $need.R) {
            $f = (& $rs -e "cat(tryCatch(as.character(packageVersion('$p')), error=function(e) ''))" 2>$null | Out-String).Trim()
            if (-not $f) { [void]$missing.Add($p) }
            $rows += [pscustomobject]@{ Item = "R: $p"; Status = $(if ($f) { 'OK' } else { 'MISSING' }); Found = $f; Tested = '' }
        }
    }
}
foreach ($t in (@('pandoc') + $need.Quick + $need.Rebuild | Select-Object -Unique)) {
    $names = if ($t -eq 'plink1.9') { @('plink1.9', 'plink') } else { @($t) }
    $c = $names | ForEach-Object { Get-Command $_ -ErrorAction SilentlyContinue } | Select-Object -First 1
    $kind = if ($t -in $need.Quick) { 'needed' } elseif ($t -eq 'pandoc') { 'optional, for -Docx' } else { 'rebuild only' }
    if (-not $c -and $kind -eq 'needed') { [void]$missing.Add($t) }
    $rows += [pscustomobject]@{ Item = "tool: $t ($kind)"; Status = $(if ($c) { 'OK' } elseif ($kind -eq 'needed') { 'MISSING' } else { 'not installed' }); Found = $(if ($c) { $c.Source } else { '' }); Tested = '' }
}
$rows | Format-Table -AutoSize | Out-String -Width 200 | Write-Host
Write-Host "Tested on: $($spec.tested_on)"
$bad = 0
foreach ($c in $Chapters) {
    $e = $spec.chapters.PSObject.Properties[[string]$c].Value
    $reqs = @($e.python) + @($spec.base.python) + @($e.r) + @($e.quick_tools)
    $lack = @($reqs | Where-Object { $missing.Contains($_) } | Select-Object -Unique)
    if ($lack.Count) { Write-Host ("Chapter {0,2}: NOT READY, missing {1}" -f $c, ($lack -join ', ')) -ForegroundColor Red; $bad++ }
    else { Write-Host ("Chapter {0,2}: ready for the quick check ({1})" -f $c, $e.title) -ForegroundColor Green }
}
exit $(if ($bad) { 1 } else { 0 })
