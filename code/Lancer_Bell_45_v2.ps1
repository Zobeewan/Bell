param(
    [double]$Coherence = 1.0,
    [int]$Paires = 200000,
    [int]$Graine = 20260916,
    [switch]$SansFenetre,
    [switch]$Rapide,
    [string]$DossierSortie = ''
)
$ErrorActionPreference = 'Stop'
$bellRoot = Split-Path -Parent $PSScriptRoot
$bellRuntime = Join-Path $env:USERPROFILE '.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$bellDeps = Join-Path $bellRoot 'tmp/mpl_runtime'
$previousPythonPath = $env:PYTHONPATH
$previousMplConfig = $env:MPLCONFIGDIR
try {
    if ((Test-Path -LiteralPath $bellRuntime) -and (Test-Path -LiteralPath $bellDeps)) {
        $bellPython = $bellRuntime
        $env:PYTHONPATH = $bellDeps + [IO.Path]::PathSeparator + $previousPythonPath
        $env:MPLCONFIGDIR = Join-Path $bellRoot 'tmp/matplotlib'
    } else {
        $bellPython = (Get-Command python -ErrorAction Stop).Source
    }
    $bellArgs = @((Join-Path $PSScriptRoot 'Bell_test_45.py'), '--coherence', $Coherence.ToString([Globalization.CultureInfo]::InvariantCulture), '--pairs', $Paires, '--seed', $Graine)
    if ($SansFenetre) { $bellArgs += '--no-show' }
    if ($Rapide) { $bellArgs += '--quick' }
    if ($DossierSortie) { $bellArgs += @('--output-dir', $DossierSortie) }
    & $bellPython @bellArgs
    $bellExitCode = $LASTEXITCODE
} finally {
    $env:PYTHONPATH = $previousPythonPath
    $env:MPLCONFIGDIR = $previousMplConfig
}
exit $bellExitCode
