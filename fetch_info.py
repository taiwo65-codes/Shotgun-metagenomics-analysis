import urllib.request
import json
import time
import re

accs = [
    ('SRX20110651', 'CRC'),
    ('SRX20110652', 'CRC'),
    ('SRX20110653', 'Healthy'),
    ('SRX20110654', 'Healthy')
]

for acc, group in accs:
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={acc}&retmode=json"
    req = urllib.request.Request(search_url, headers={"User-Agent": "BioPython/1.0"})
    with urllib.request.urlopen(req) as resp:
        sdata = json.loads(resp.read().decode())
    uid = sdata["esearchresult"]["idlist"][0]
    time.sleep(0.5)
    sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id={uid}&retmode=json"
    req2 = urllib.request.Request(sum_url, headers={"User-Agent": "BioPython/1.0"})
    with urllib.request.urlopen(req2) as resp:
        data = json.loads(resp.read().decode())
    item = data["result"][uid]
    runs = item.get("runs", "")
    match = re.search(r'acc="([A-Z0-9]+)"[^>]*total_spots="([0-9]+)"[^>]*total_bases="([0-9]+)"', runs)
    if match:
        run_acc, spots, bases = match.groups()
        gb = int(bases) / 1e9
        print(f"Accession: {acc} | Group: {group} | Run: {run_acc} | Spots: {spots} | Bases: {bases} (~{gb:.2f} Gbp)")
    else:
        print(f"Accession: {acc} | Group: {group} | Runs raw: {runs}")
    time.sleep(0.5)
