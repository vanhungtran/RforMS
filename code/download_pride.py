#!/usr/bin/env python3
"""List and download PRIDE datasets for R for Mass Spectrometry book chapter."""
import urllib.request
import json
import os
import sys
import subprocess

DATASETS = {
    "PXD004886": "DIA benchmark (Bruderer et al.) - known ground truth",
    "PXD010154": "Plasma proteomics, multiple search modes",
    "PXD000547": "DIA/SWATH paired clinical samples",
}

BASE = "/cluster/data/ck_care/ltran/r4ms_book/raw"
PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v1"


def pride_request(endpoint):
    """Make a request to PRIDE API with JSON accept header."""
    url = f"{PRIDE_API}/{endpoint}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = resp.read().decode()
            print(f"  DEBUG: response length={len(raw)}, status={resp.status}")
            if len(raw) == 0:
                print(f"  DEBUG: Empty response body")
                return None
            return json.loads(raw)
    except Exception as e:
        print(f"  API error for {url}: {e}")
        return None


def list_files(accession):
    """List files for a PRIDE project."""
    data = pride_request(f"files/byProject?accession={accession}")
    if not data:
        return []
    return data


def download_file(url, dest_dir, fname):
    """Download a single file using wget."""
    dest = os.path.join(dest_dir, fname)
    if os.path.exists(dest):
        print(f"    Already exists: {fname}")
        return True
    print(f"    Downloading: {fname}")
    result = subprocess.run(
        ["wget", "-q", "--show-progress", "-O", dest, url],
        capture_output=True, text=True, timeout=300
    )
    return result.returncode == 0


def get_project_info(acc):
    """Get project metadata including submission date."""
    data = pride_request(f"projects/{acc}")
    if not data:
        return None
    sub_date = data.get("submissionDate", "?")
    pub_date = data.get("publicationDate", "?")
    refs = data.get("references", [])
    pmid = refs[0].get("pubmedId", "?") if refs else "?"
    return {"submissionDate": sub_date, "publicationDate": pub_date, "pmid": pmid}


def ftp_download(acc, dest_dir):
    """Try to download from PRIDE FTP using known year/month from project date."""
    info = get_project_info(acc)
    if not info:
        return False
    
    sub_date = info["submissionDate"]
    if sub_date and sub_date != "?":
        year = sub_date[:4]
        # Try several month patterns
        months_to_try = [sub_date[5:7]]
        # Also try adjacent months
        m = int(sub_date[5:7])
        months_to_try.extend([f"{x:02d}" for x in [m-1, m+1] if 1 <= x <= 12])
        
        for month in months_to_try:
            ftp_url = f"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/{year}/{month}/{acc}/"
            print(f"  Trying FTP: {ftp_url}")
            result = subprocess.run(
                ["wget", "--no-passive-ftp", "-q", "--spider", ftp_url],
                capture_output=True, text=True, timeout=15
            )
            if result.returncode == 0:
                print(f"  Found at: {ftp_url}")
                # Download non-raw files
                subprocess.run([
                    "wget", "--no-passive-ftp", "-r", "-np", "-nH", "--cut-dirs=5",
                    "-A", "*.txt,*.csv,*.tsv,*.xlsx,*.pdf,*.fasta,*.xml,*.mzid,*.mzTab",
                    "-P", dest_dir, ftp_url
                ])
                return True
    
    print(f"  Could not find FTP path for {acc}")
    return False


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "list"

    for acc, desc in DATASETS.items():
        print(f"\n{'='*60}")
        print(f"  {acc}: {desc}")
        print(f"{'='*60}")

        # Show project metadata
        info = get_project_info(acc)
        if info:
            print(f"  Submitted: {info['submissionDate']}  |  PMID: {info['pmid']}")

        # Try PRIDE REST API for files
        files = list_files(acc)
        
        if not files:
            print("  REST API returned no files. Trying FTP...")
            dest_dir = os.path.join(BASE, acc)
            os.makedirs(dest_dir, exist_ok=True)
            if mode == "download":
                ftp_download(acc, dest_dir)
            else:
                # Just show what we'd try
                if info and info["submissionDate"] and info["submissionDate"] != "?":
                    y, m = info["submissionDate"][:4], info["submissionDate"][5:7]
                    print(f"  Would try: ftp://ftp.pride.ebi.ac.uk/pride/data/archive/{y}/{m}/{acc}/")
            continue

        total_mb = sum(f.get("fileSize", 0) for f in files) / (1024 * 1024)
        print(f"  {len(files)} files, {total_mb:.0f} MB total\n")

        for f in files:
            name = f.get("fileName", "?")
            size = f.get("fileSize", 0) / (1024 * 1024)
            link = f.get("downloadLink", "")
            print(f"  {size:8.1f} MB  {name}")

        if mode == "download":
            dest_dir = os.path.join(BASE, acc)
            os.makedirs(dest_dir, exist_ok=True)

            # Only download processed/metadata files, skip raw instrument data
            skip_ext = (".raw", ".wiff", ".wiff.scan", ".d", ".mgf")
            for f in files:
                name = f.get("fileName", "")
                link = f.get("downloadLink", "")
                if not name or not link:
                    continue
                if name.endswith(skip_ext):
                    print(f"    Skipping raw: {name}")
                    continue
                download_file(link, dest_dir, name)

            print(f"\n  Downloaded to: {dest_dir}")


if __name__ == "__main__":
    main()
