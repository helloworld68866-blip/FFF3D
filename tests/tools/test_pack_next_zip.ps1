Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Assert-True {
    param(
        [bool]$Condition,
        [string]$Message
    )

    if (-not $Condition) {
        throw $Message
    }
}

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$scriptPath = Join-Path $repoRoot 'pack-next-zip.ps1'
$exePath = Join-Path $repoRoot 'pack-next-zip.exe'
$analysisExePath = Join-Path $repoRoot 'pack-next-zip-analysis.exe'

$tempRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("dec3d-pack-test-{0}" -f [System.Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $tempRoot | Out-Null

try {
    foreach ($directoryName in @('artifacts', 'cases', 'docs', 'src', 'tests')) {
        $directoryPath = Join-Path $tempRoot $directoryName
        New-Item -ItemType Directory -Path $directoryPath | Out-Null
        Set-Content -LiteralPath (Join-Path $directoryPath 'sample.txt') -Value $directoryName
    }

    $analysisPath = Join-Path $tempRoot 'analysis'
    New-Item -ItemType Directory -Path $analysisPath | Out-Null
    Set-Content -LiteralPath (Join-Path $analysisPath 'sample.txt') -Value 'analysis'

    Set-Content -LiteralPath (Join-Path $tempRoot 'AGENTS.md') -Value 'agents'
    Set-Content -LiteralPath (Join-Path $tempRoot 'README.md') -Value 'readme'
    Set-Content -LiteralPath (Join-Path $tempRoot 'CMakeLists.txt') -Value 'cmake'
    Set-Content -LiteralPath (Join-Path $tempRoot '7.zip') -Value 'old-zip'
    Set-Content -LiteralPath (Join-Path $tempRoot '40.zip') -Value 'old-zip'
    Set-Content -LiteralPath (Join-Path $tempRoot 'notes.zip') -Value 'ignored'

    $output = & powershell -NoProfile -ExecutionPolicy Bypass -File $scriptPath -RootPath $tempRoot 2>&1
    $exitCode = $LASTEXITCODE

    Assert-True ($exitCode -eq 0) ("Expected pack script to succeed, but it failed.`n{0}" -f ($output -join [Environment]::NewLine))

    $archivePath = Join-Path $tempRoot '41.zip'
    Assert-True (Test-Path -LiteralPath $archivePath) "Expected 41.zip to be created."

    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $archive = [System.IO.Compression.ZipFile]::OpenRead($archivePath)
    try {
        $entryNames = $archive.Entries | ForEach-Object { $_.FullName.Replace('\', '/') }

        foreach ($expectedEntry in @(
            'artifacts/sample.txt',
            'cases/sample.txt',
            'docs/sample.txt',
            'src/sample.txt',
            'tests/sample.txt',
            'AGENTS.md',
            'README.md',
            'CMakeLists.txt'
        )) {
            Assert-True ($entryNames -contains $expectedEntry) ("Expected archive to include '{0}', but entries were:`n{1}" -f $expectedEntry, ($entryNames -join [Environment]::NewLine))
        }

        Assert-True (-not ($entryNames -contains 'pack-next-zip.ps1')) "The archive should not contain the pack script itself."
        Assert-True (-not ($entryNames -contains 'analysis/sample.txt')) "The archive should not contain the analysis directory."
    }
    finally {
        $archive.Dispose()
    }

    Assert-True (Test-Path -LiteralPath $exePath) "Expected the no-window launcher pack-next-zip.exe to exist."
    Assert-True (Test-Path -LiteralPath $analysisExePath) "Expected the analysis launcher pack-next-zip-analysis.exe to exist."

    Copy-Item -LiteralPath $exePath -Destination (Join-Path $tempRoot 'pack-next-zip.exe')
    Copy-Item -LiteralPath $analysisExePath -Destination (Join-Path $tempRoot 'pack-next-zip-analysis.exe')

    $launcherOutput = & (Join-Path $tempRoot 'pack-next-zip.exe') 2>&1
    $launcherExitCode = $LASTEXITCODE

    Assert-True ($launcherExitCode -eq 0) ("Expected no-window launcher to succeed, but it failed.`n{0}" -f ($launcherOutput -join [Environment]::NewLine))
    Assert-True (Test-Path -LiteralPath (Join-Path $tempRoot '42.zip')) "Expected no-window launcher to create 42.zip."

    $analysisLauncherOutput = & (Join-Path $tempRoot 'pack-next-zip-analysis.exe') 2>&1
    $analysisLauncherExitCode = $LASTEXITCODE

    Assert-True ($analysisLauncherExitCode -eq 0) ("Expected analysis launcher to succeed, but it failed.`n{0}" -f ($analysisLauncherOutput -join [Environment]::NewLine))

    $analysisArchivePath = Join-Path $tempRoot '43.zip'
    Assert-True (Test-Path -LiteralPath $analysisArchivePath) "Expected analysis launcher to create 43.zip."

    $analysisArchive = [System.IO.Compression.ZipFile]::OpenRead($analysisArchivePath)
    try {
        $analysisEntryNames = $analysisArchive.Entries | ForEach-Object { $_.FullName.Replace('\', '/') }
        Assert-True ($analysisEntryNames -contains 'analysis/sample.txt') "Expected analysis archive to include 'analysis/sample.txt'."
    }
    finally {
        $analysisArchive.Dispose()
    }

    Write-Host 'PASS: pack-next-zip creates the next numbered archive with the requested contents.'
}
finally {
    Remove-Item -LiteralPath $tempRoot -Recurse -Force -ErrorAction SilentlyContinue
}
