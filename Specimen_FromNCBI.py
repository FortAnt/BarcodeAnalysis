#! /usr/bin/env python3.7

#usage python Specimen_FromNCBI.py <NCBI_accession_numbers_file.txt> > NCBI_specimens.txt
from Bio import Entrez
import sys
inFile = sys.argv[1]



accessions_file = inFile
with open(accessions_file) as f:
    ids = f.read().split('\n')

# Fetch the entries from Entrez
Entrez.email = ''  # Insert your NCBI email here
handle = Entrez.efetch('nuccore', id=ids, retmode='xml')
response = Entrez.read(handle)

def extract_strain(entry):
    sources = [feature for feature in entry['GBSeq_feature-table']
               if feature['GBFeature_key'] == 'source']

    for source in sources:
        qualifiers = [qual for qual in source['GBFeature_quals']
                      if qual['GBQualifier_name'] == 'specimen_voucher' or qual['GBQualifier_name'] == 'strain' or qual['GBQualifier_name'] == 'isolate']

        for qualifier in qualifiers:
            yield qualifier['GBQualifier_value']

for entry in response:
    accession = entry['GBSeq_primary-accession']
    for strain in extract_strain(entry):
            print(accession,strain, sep='\t')
