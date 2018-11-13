# clinvar_snps
Annotating pathogenic ClinVar SNPs with nearby PAM sequences.

Author: Ruth Hanna
Email: rhanna@broadinstitute.org
Title: Identifying PAMs near ClinVar SNPs

Description: This script takes in 3 inputs:
1. A raw file of SNP positions downloaded from the ClinVar website (variant_summary.txt):
ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
This file is updated weekly.

2. A PAM file with a list of Cas nucleases and information about their PAMs and editing windows.
Two sample PAM files, cas_nuclease_pams.txt and base_editor_pams.txt are provided here. These can be used to generate the number of editable SNPs with all current Cas nucleases and base editors, respectively.

3. The desired name of the output file.

To run, execute the following in the terminal:

python annnotate_clinvar_snps.py --variant_file variant_summary.txt --pam_file [pam_file.txt] --output_name [outputname]

where [pam_file.txt] is the name of the PAM file and [outputname] is the desired folder name (e.g. base_editor_snps).
