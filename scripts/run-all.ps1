<#
.SYNOPSIS  Run the quick check (or the full rerun) of several chapters and print one summary table.
.EXAMPLE   pwsh scripts/run-all.ps1 -Chapters 1,2,4,12,13
#>
[CmdletBinding()]
param([string[]]$Chapters = @('1-13'), [switch]$Full)
. (Join-Path $PSScriptRoot '_common.ps1')
$Chapters = ConvertTo-ChapterList $Chapters
$res = @()
foreach ($c in $Chapters) {
    $sw = [Diagnostics.Stopwatch]::StartNew()
    $p = @{ Chapter = $c }; if ($Full) { $p['Full'] = $true }
    & (Join-Path $PSScriptRoot 'run-chapter.ps1') @p | Out-Host
    $code = $LASTEXITCODE; $sw.Stop()
    $res += [pscustomobject]@{ Chapter = $c; Result = $(if ($code -eq 0) { 'PASS' } else { 'FAIL' }); Seconds = [int]$sw.Elapsed.TotalSeconds }
}
$res | Format-Table -AutoSize | Out-String | Write-Host
exit $(if ($res | Where-Object Result -eq 'FAIL') { 1 } else { 0 })
