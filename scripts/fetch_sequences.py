#!/usr/bin/env python3
import sys
import requests
import time

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def fetch_all(database, accessions_file, output_file="downloaded.fasta"):
    accessions = open(accessions_file).read().splitlines()
    accessions = [a.strip() for a in accessions if a.strip()]
    print(f"Found {len(accessions)} accession numbers")

    with open(output_file, "w") as out:
        for i, acc in enumerate(accessions, 1):
            print(f"  [{i}/{len(accessions)}] {acc} ...", end=" ", flush=True)
            try:
                r = requests.get(
                    BASE_URL,
                    params={
                        "db": database,
                        "id": acc,
                        "rettype": "fasta",
                        "retmode": "text"
                    },
                    timeout=30
                )
                r.raise_for_status()
                if r.text.startswith(">"):
                    out.write(r.text.strip() + "\n\n")
                    print("OK")
                else:
                    print("WARNING: unexpected response")
            except Exception as e:
                print(f"ERROR: {e}")
            time.sleep(0.4)

    print(f"\nSaved to: {output_file}")
    return output_file

# Alias for pipeline compatibility
def fetch_fasta(database, accessions_file, output_file="downloaded.fasta"):
    return fetch_all(database, accessions_file, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fetch_sequences.py <database> <file.txt>")
        sys.exit(1)
    fetch_all(sys.argv[1], sys.argv[2])

