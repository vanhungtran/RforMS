# Download MetaboLights mzML files for RforMS book chapters
# Handles FTP filename encoding for special characters (commas, parentheses, etc.)
# Usage: .\download_mtbls.ps1 [-DryRun] [-Dataset MTBLS38,MTBLS1455,MTBLS234]
param(
    [switch]$DryRun,
    [string[]]$Datasets = @("MTBLS38","MTBLS1455","MTBLS234")
)

$ErrorActionPreference = "Continue"
$FTP_BASE = "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public"
$LOCAL_BASE = "C:\server_Cardiocare\r4ms_book\raw"

# URL-encode filenames with special characters for FTP
function Format-FtpUrl {
    param([string]$BasePath, [string]$FileName)
    $encoded = [uri]::EscapeDataString($FileName)
    # Fix over-encoding: some chars should stay literal for FTP
    return "$BasePath/$encoded"
}

foreach ($acc in $Datasets) {
    Write-Host "`n=== $acc ===" -ForegroundColor Cyan
    $dest = Join-Path $LOCAL_BASE $acc
    New-Item -ItemType Directory -Force -Path $dest | Out-Null

    $ftpDir = "$FTP_BASE/$acc/FILES/"

    # Get file listing via .NET WebClient (handles FTP better than curl on Windows)
    $request = [System.Net.FtpWebRequest]::Create($ftpDir)
    $request.Method = [System.Net.WebRequestMethods+Ftp]::ListDirectory
    $request.Timeout = 30000

    try {
        $response = $request.GetResponse()
        $stream = $response.GetResponseStream()
        $reader = New-Object System.IO.StreamReader($stream)
        $listing = $reader.ReadToEnd() -split "`r`n" | Where-Object { $_ -ne "" }
        $reader.Close()
        $response.Close()
    } catch {
        Write-Host "  FTP listing failed: $_" -ForegroundColor Red
        continue
    }

    # Filter: mzML files only, skip raw/directories
    $mzmlFiles = $listing | Where-Object { $_ -match '\.mzML$' -and $_ -notmatch '^d' }
    $metaFiles = $listing | Where-Object { $_ -match '\.(tsv|csv|txt|xlsx)$' }

    Write-Host "  Found: $($mzmlFiles.Count) mzML, $($metaFiles.Count) metadata files" -ForegroundColor Green

    # Download metadata first (always small)
    foreach ($f in $metaFiles) {
        $outFile = Join-Path $dest $f
        if ($DryRun) {
            Write-Host "  WOULD download: $f"
        } else {
            try {
                $url = Format-FtpUrl -BasePath $ftpDir -FileName $f
                (New-Object System.Net.WebClient).DownloadFile($url, $outFile)
                Write-Host "  OK: $f" -ForegroundColor Green
            } catch {
                Write-Host "  FAIL: $f" -ForegroundColor Red
            }
        }
    }

    # Download mzML files
    $count = 0
    $total = $mzmlFiles.Count
    foreach ($f in $mzmlFiles) {
        $count++
        $outFile = Join-Path $dest $f

        if (Test-Path $outFile) {
            $existingSize = (Get-Item $outFile).Length
            if ($existingSize -gt 1000) {
                Write-Host "  [$count/$total] SKIP (exists): $f"
                continue
            }
        }

        if ($DryRun) {
            Write-Host "  [$count/$total] WOULD: $f"
        } else {
            try {
                $url = Format-FtpUrl -BasePath $ftpDir -FileName $f
                Write-Host "  [$count/$total] DOWNLOAD: $f ..." -NoNewline
                (New-Object System.Net.WebClient).DownloadFile($url, $outFile)
                $size = (Get-Item $outFile).Length / 1MB
                Write-Host " $([math]::Round($size,1)) MB" -ForegroundColor Green
            } catch {
                Write-Host " FAILED: $_" -ForegroundColor Red
            }
        }
    }
}

Write-Host "`n=== Complete ===" -ForegroundColor Cyan
if ($DryRun) { Write-Host "Run without -DryRun to download." }
