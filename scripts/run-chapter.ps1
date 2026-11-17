<#
.SYNOPSIS  Check one chapter or run one exercise from the shipped data.
.EXAMPLE
  pwsh scripts/run-chapter.ps1 -Chapter 12                        # recompute every number of the chapter
  pwsh scripts/run-chapter.ps1 -Chapter 12 -Full                  # also rerun every script and compare each result table byte for byte
  pwsh scripts/run-chapter.ps1 -Chapter 13 -Script ex13_5_sir_gi.py   # run one exercise; the output stays in chapters/chapter13/.work
  pwsh scripts/run-chapter.ps1 -Chapter 6 -Full -Sim              # include the long TreeTime reruns of Chapter 6
  pwsh scripts/run-chapter.ps1 -Chapter 12 -Docx C:\path\Chapter_12.docx   # also compare your copy of the chapter text (not in the repository)
#>
[CmdletBinding()]
param([Parameter(Mandatory)][ValidateRange(1, 13)][int]$Chapter, [switch]$Full, [switch]$Doi, [switch]$Sim, [string]$Script, [string]$Docx, [switch]$Offline)
. (Join-Path $PSScriptRoot '_common.ps1')
Set-RunEnvironment
$py = Get-Python
if (-not $py) { throw 'Python was not found. Run scripts/setup.ps1 first.' }
$a = @((Join-Path $RepoRoot 'tools/run_chapter.py'), $Chapter)
if ($Full) { $a += '--full' }; if ($Doi) { $a += '--doi' }; if ($Sim) { $a += '--sim' }; if ($Offline) { $a += '--offline' }
if ($Script) { $a += @('--script', $Script) }
if ($Docx) { $a += @('--docx', (Resolve-Path $Docx).Path) }
Push-Location $RepoRoot
try { & $py @a; exit $LASTEXITCODE } finally { Pop-Location }
