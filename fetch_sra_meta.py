import urllib.request
import xml.etree.ElementTree as ET

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRX20110651,SRX20110652,SRX20110653,SRX20110654&retmode=xml'
req = urllib.request.Request(url, headers={'User-Agent': 'BioPython/1.0'})
with urllib.request.urlopen(req) as resp:
    xml_data = resp.read()

root = ET.fromstring(xml_data)
for exp in root.findall('.//EXPERIMENT_PACKAGE'):
    acc = exp.find('.//EXPERIMENT').attrib.get('accession')
    title_elem = exp.find('.//EXPERIMENT/TITLE')
    title = title_elem.text if title_elem is not None else ''
    bp_elem = exp.find('.//EXTERNAL_ID[@namespace="BioProject"]')
    bioproject = bp_elem.text if bp_elem is not None else ''
    bs_elem = exp.find('.//EXTERNAL_ID[@namespace="BioSample"]')
    biosample = bs_elem.text if bs_elem is not None else ''
    run_elem = exp.find('.//RUN')
    run = run_elem.attrib.get('accession') if run_elem is not None else ''
    org_elem = exp.find('.//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME')
    organism = org_elem.text if org_elem is not None else ''
    
    # Extract sample attributes
    attrs = {}
    for attr in exp.findall('.//SAMPLE_ATTRIBUTE'):
        tag = attr.find('TAG')
        val = attr.find('VALUE')
        if tag is not None and val is not None and tag.text and val.text:
            attrs[tag.text] = val.text

    print(f"Exp: {acc} | Run: {run} | BioProject: {bioproject} | BioSample: {biosample} | Organism: {organism}")
    print(f"  Title: {title}")
    print(f"  Attributes: {attrs}")
