# BarcodeAnalysis

This repository contains the scripts and command lines used to generate phylogenetic trees and GMYC species delimitation analysis.

This pipeline was tested on Ubuntu 20.04, R version 3.5.1 and Python version 3.7.
The python scripts need to be manually modified to add your personal NCBI registered email address.

**1) Preparation of alignments**
======
Download sequences from NCBI, to obtain "barcode.fasta"
For example using *Ulva [organism] AND rbcL [gene] AND plastid [filter]* in the NCBI (https://www.ncbi.nlm.nih.gov/nuccore/),
then click "send to" > "file" > "fasta" > "create file". One can apply a length filter to omit full organellar genomes.

**1.1) For alignments retaining all fasta headers from NCBI, and removing un-characterised sequences containing Genus "sp." :**
------
This keeps all headers (including description of the sequence).\
Replace spaces by "_" from voucher names, required for downstream analysis:

>sed -i "s/ /_/g" barcode.fasta

Obtain voucher names that do not contain "Genus sp":

>grep '>' barcode.fasta | grep -v 'Genus_sp.' > NO_Genus_sp.IDs

Remove ">" from the list, required for samtools faidx

>sed -i "s/>//g" NO_Genus_sp.IDs

Extract sequences matching the NO_Genus_sp.IDs list using samtools, available here: http://www.htslib.org/

>xargs samtools faidx barcode.fasta < NO_Genus_sp.IDs > barcode_noGenussp.fasta

**1.2) For Alignments with fasta headers in the format: NCBIaccessionnumber_Genus_species, retaining all species names, even uncharacterized ones**
------
This leads to much cleaner fasta headers, retaining only Accession number and species information, delimited by "_".\
Obtain accession numbers:

>grep '>\S*' barcode.fasta | awk '{print $1}' | sed -e "s/>//g" > NCBI_vouchers.txt
 
Run python script to obtain "organism" metadata for each accessions, then sort accessions alphabetically:\
It requires the entrez python API, available here: https://biopython.org/wiki/Download to retrieve metadata from the accession numbers.
>python [Species_FromNCBI.py](https://github.com/FortAnt/BarcodeAnalysis/blob/main/Species_FromNCBI.py) NCBI_vouchers.txt | sort > species_NCBI.txt

Reformat python hits by concatenating accession number and species name, separated by "_":

>sed -e "s|\t|\_|g" species_NCBI.txt | sed -e "s/ \/_/g" > renamed_vouchers.txt

Get full fasta headers, sort alphabetically:

>grep '>' barcode.fasta | sed -e "s/>//g" | sort > full_fasta_headers.txt

Concatenate fasta headers and renamed headers file

>paste -d"\t" full_fasta_headers.txt renamed_vouchers.txt > renamed_headers.txt

Rename fasta headers with new headers, using seqkit (https://github.com/shenwei356/seqkit):

>seqkit replace -p "(.+)" -r '{kv}' -k ./renamed_headers.txt ./barcode.fasta > barcode.renamed.fasta

**1.3) Alignment**
------
Align sequences, using mafft available here: https://mafft.cbrc.jp/alignment/software/

The --adjustdirection helps as some NCBI sequences are reverse complement compared to others.\
The --maxiterate 1000 is optional but is used for alignments containing gaps or between very divergent sequences.

>mafft --adjustdirection --maxiterate 1000 --threads X barcode.renamed.fasta > barcode.renamed.aln.fasta

Remove positions with more than X unknowns, using trimal available here: http://trimal.cgenomics.org/downloads

>trimal -in barcode.renamed.aln.fasta -out barcode.renamed.trimmed_gtX.aln.fasta -gt X

Remove sequences with less than Y overlap:

>trimal -in barcode.renamed.trimmed_gtX.aln.fasta -out barcode.renamed.trimmed_gtXseqY.aln.fasta -resoverlap 0.70 -seqoverlap Y 

If no gaps expected between sequences (such as COI or rbcL barcode), replace gaps with Ns:

>sed -i "s/-/n/g" barcode.renamed.trimmed_gtXseqY.aln.fasta

Alternatively, if gaps are present within the sequence, manually replace gaps in 5’ and 3’ end in the alignment directly in Jalview, available here: https://www.jalview.org/

**2) Phylogenetic analysis**
======
•	Use jModeltest2 for selecting best evolutionary model. Download here: https://github.com/ddarriba/jmodeltest2/releases

•	Select best evolutionary model based on AIB and/or BIC score.

•	Maximum likelihood analysis using raxml-ng, available here: https://github.com/amkozlov/raxml-ng.

•	Bayesian probability tree using mrbayes (MPI version for faster analysis), available here: https://nbisweden.github.io/MrBayes/download.html.

•	GMYC analysis using BEAUti/BEAST and the BEAGLE library for faster computation, available here: https://beast.community/programs. We found this tutorial useful: https://francoismichonneau.net/gmyc-tutorial/

•	Summarise the BEAST trees with treannotator & suitable burnin

> treeannotator -burnin X BEAST.trees > BEAST.trees.nex 

•	Species delimitation and obtain the percentage of agreement in species names using R.

>[R script available here](https://github.com/FortAnt/BarcodeAnalysis/blob/main/GMYC_delimitation_script.R)

•	Visualise trees using Figtree, available here: http://tree.bio.ed.ac.uk/software/figtree/.

**3) Distribution analysis**
=====
•	Obtain the list of accession numbers from the alignment, removing ">" and retaining only accession numbers

>grep '>'  barcode.renamed.trimmed_gtXseqY.aln.fasta | sed "s/\..*//g" | sed "s/>//g" > accession_numbers_barcodes.txt

•	Use the entrez python API, available here: https://biopython.org/wiki/Download to retrieve metadata from the accession numbers. Scripts are included in this repository  and were adapted from Juan Manuel Berros contribution to https://www.biostars.org/p/221662/.

>python [GPS_coordinates_FromNCBI.py](https://github.com/FortAnt/BarcodeAnalysis/blob/main/Specimen_FromNCBI.py) accession_numbers_barcodes.txt > NCBI_GPS.txt

>python [Specimen_FromNCBI.py](https://github.com/FortAnt/BarcodeAnalysis/blob/main/Specimen_FromNCBI.py) accession_numbers_barcodes.txt > NCBI_specimens.txt

•	Use R (or excel) for merging GPS, specimens, species and remove duplicate specimens. Then select species of interest and plot distribution using the Rworldmap package in R (https://www.rdocumentation.org/packages/rworldmap/versions/1.3-6)



**How to modify metadata on your own NCBI records**
=====
When nomenclature changes for your species of interest, of if one's own record was misannotated, it is important to modify your records with the latest appropriate nomenclature. This avoids the amplification of misannotations in the database.

Is is incredibly simple to [update the metadata associated with a NCBI record](https://www.ncbi.nlm.nih.gov/genbank/update/#:~:text=All%20update%20files%20should%20be,to%20update%20an%20existing%20record)

To do so, from the email account associated with your NCBI submissions, send an email to gb-admin@ncbi.nlm.nih.gov, with as attachment a **tab-delimited** text file containing in the first column the accession numbers (the column is named acc. num.), and the other columns containing the metadata to correct.

For example, to modify the "organism" information of some of my vouchers that were wrongly annotated as U. laetevirens, the text file contained:


acc. num. organism\
MT160559  Ulva lacinulata\
MT160599  Ulva lacinulata\
MT160596  Ulva lacinulata\
MT160544  Ulva lacinulata\
MT160614  Ulva lacinulata\
MT160549  Ulva lacinulata\
MT160548  Ulva lacinulata\
MT160551  Ulva lacinulata


To modify additional metadata information such as strain name, cultivar, GPS coordinates, one can add as many columns to the text file as metadata to modify. [Here](https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/) is the list of source modifiers of the NCBI.
