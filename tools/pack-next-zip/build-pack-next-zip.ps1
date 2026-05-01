Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$scriptDirectory = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent (Split-Path -Parent $scriptDirectory)
$sourcePath = Join-Path $scriptDirectory 'Program.cs'
$sourceCode = Get-Content -Raw -LiteralPath $sourcePath

function Build-PackNextZipExe {
    param(
        [string]$OutputFileName,
        [string]$CompiledSourceCode
    )

    $outputPath = Join-Path $repoRoot $OutputFileName

    if (Test-Path -LiteralPath $outputPath) {
        Remove-Item -LiteralPath $outputPath -Force
    }

    Add-Type `
        -Language CSharp `
        -TypeDefinition $CompiledSourceCode `
        -OutputAssembly $outputPath `
        -OutputType WindowsApplication `
        -ReferencedAssemblies @(
            'System.IO.Compression',
            'System.IO.Compression.FileSystem',
            'System.Windows.Forms'
        )
}

$defaultSourceCode = $sourceCode.Replace('/*ANALYSIS_DIRECTORY*/', '')
$analysisSourceCode = $sourceCode.Replace('/*ANALYSIS_DIRECTORY*/', ',' + [Environment]::NewLine + '        "analysis"')

Build-PackNextZipExe -OutputFileName 'pack-next-zip.exe' -CompiledSourceCode $defaultSourceCode
Build-PackNextZipExe -OutputFileName 'pack-next-zip-analysis.exe' -CompiledSourceCode $analysisSourceCode
