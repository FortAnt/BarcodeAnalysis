#! /usr/bin/env python3.7


#usage python GPS_coordinates_FromNCBI.py <NCBI_accession_numbers_file.txt> > NCBI_GPS.txt
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


def extract_latlon(entry):
    sources = [feature for feature in entry['GBSeq_feature-table']
               if feature['GBFeature_key'] == 'source']

    for source in sources:
        qualifiers = [qual for qual in source['GBFeature_quals']
                      if qual['GBQualifier_name'] == 'lat_lon']

        for qualifier in qualifiers:
            yield qualifier['GBQualifier_value']


for entry in response:
    accession = entry['GBSeq_primary-accession']
    for lat_lon in extract_latlon(entry):
            print(accession, lat_lon, sep='\t')