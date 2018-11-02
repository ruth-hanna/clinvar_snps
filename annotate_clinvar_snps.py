"""
Author: Ruth Hanna
Email: rhanna@broadinstitute.org
Title: Identifying PAMs near ClinVar SNPs
Description: This script reads in two input files:
1. A raw file of SNP positions downloaded from the ClinVar website.
2. A file with a list of Cas nucleases and information about their PAMs and editing windows.
The script first identifies the genomic context in the reference genome for each SNP and then searches for PAMs for each Cas in the genomic context.
The output is a list of sgRNAs, annotated with the SNP and Cas, and a list of SNPs, annotated with whether they are editable by each Cas.
All genomic locations use 1-based indexing; all locations in sequences (e.g. PAM_Location) use 0-based indexing.
"""

import pandas as pd
import numpy as np
import argparse, csv, itertools, gzip, os
from datetime import datetime

def GetParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variant_file',
        type=str,
        help='.txt input file from the ClinVar website)')
    parser.add_argument('--pam_file',
    	type=str,
    	help='.txt input file with the following columns: name of Cas, PAM, start of editing window, end of editing window, sgRNA length, allowed nt (e.g. for base editing)')
    parser.add_argument('--output_name',
    	type =str,
    	help='name of output file')
    return parser

def ParseVariantFile(variant_df):
	# Cleans variant file from ClinVar database.
	# Removes non-GRCh38 rows, non-SNPs, and chromosomes other than 1-22, X, and Y.
	# Returns cleaned df.
	
	print('Parsing variant file')
	temp_variant_df = variant_df[(variant_df.Assembly == 'GRCh38') & (variant_df.ClinicalSignificance == 'Pathogenic') & (variant_df.Type == 'single nucleotide variant') & (variant_df.Chromosome != 'MT') & (variant_df.Chromosome != 'na') & (variant_df.ReferenceAllele != 'na')]
	parsed_variant_df = temp_variant_df.copy()
	parsed_variant_df.index = range(0,len(parsed_variant_df))
	parsed_variant_df.rename(columns = {'Start':'ClinVar_SNP_Position'}, inplace=True)
	return parsed_variant_df

def ParsePAMFile(pam_df):
	# Converts input PAMs from a string to a list.
	split_row_list = []
	column_list = list(pam_df.columns.values)
	column_list[column_list.index('PAM')] = 'Split_PAM'
	for row in pam_df['PAM']:
		split_row = row.split(',')
		split_row_list.append(split_row)
	pam_df['Split_PAM'] = split_row_list
	pam_df = pam_df.drop(['PAM'],axis=1)
	pam_df = pam_df[column_list]
	return pam_df

def ReadChromFiles(start_chrom,end_chrom):
	# Loads raw chromosome files and stores as a dictionary.

	chr_dict = {}
	for chrom_num in range(start_chrom,end_chrom+1):
		if chrom_num == 23:
			chrom_num = 'X'
		if chrom_num == 24:
			chrom_num = 'Y'
		fasta_path = '/rnai/refdb/genome/human/GRCh38/'
		chr_file = fasta_path + 'hs_ref_GRCh38_chr' + str(chrom_num) + '.fa.gz'
		print('Reading chrom ' + str(chrom_num))
		f = gzip.open(chr_file,'r')
		chr_seq = f.read()
		f.close()
		chr_seq_split = chr_seq.split('\n')
		chr_seq = ''.join(chr_seq_split[1:])
		chr_dict[str(chrom_num)] = chr_seq
	return chr_dict

def LocateSNP(chr_dict, parsed_variant_df):
	# Takes in a df of SNP genomic locations and the dictionary of chromosome sequences.
	# For each SNP, extracts the 50 nts preceding the SNP and the 50 nts following the SNP.
	# Outputs a snp_df with the genomic context of the given SNP.

	snp_location_list = []
	error_list = []
	snp_df_columns_list = list(parsed_variant_df.columns.values)
	snp_df_columns_list.extend(['Genomic_Context_Start','Genomic_Context_End','SNP_Position','Preceding_Genomic_Context','Succeeding_Genomic_Context','Genomic_Context'])
	for index, row in parsed_variant_df.iterrows():
		print('Extracting genomic context for row ' + str(index))
		r = list(row)
		snp_chr_seq = chr_dict[row['Chromosome']]
		snp_position = row['ClinVar_SNP_Position']
		# Extract surrounding genomic context (taking into account 0-based indexing in Python).
		preceding_genomic_context = snp_chr_seq[(snp_position-51):snp_position-1]
		succeeding_genomic_context = snp_chr_seq[(snp_position):(snp_position+50)]
		snp_genomic_context = preceding_genomic_context + row['AlternateAllele'] + succeeding_genomic_context
		r.extend([snp_position-50, snp_position+50, snp_position, preceding_genomic_context, succeeding_genomic_context, snp_genomic_context])
		snp_location_list.append(r)
		# Check that reference alleles match. If they don't, append to error_list.
		reference_allele_from_chrom_seq = snp_chr_seq[snp_position-1]
		if reference_allele_from_chrom_seq != row['ReferenceAllele']:
			print('Reference alleles do not match for SNP # ' + row['#AlleleID'])
			error_list.append(r)
		# Check that length of extracted genomic sequence is correct. If it isn't, append to error_list.
		if len(snp_genomic_context) != 101:
			print('Error: genomic context is not long enough for snp = ' + row['#AlleleID'])
			error_list.append(r)
	snp_df = pd.DataFrame(data=snp_location_list,columns = snp_df_columns_list)
	if error_list:
		snp_error_df = pd.DataFrame(data=error_list,columns = snp_df_columns_list)
	if not error_list:
		snp_error_df = pd.DataFrame(data=[],columns=[])
	return (snp_error_df, snp_df)
	
def ReverseComplement(sequence):
	# Returns the reverse complement of input sequence.

	complements_dict = {'A':'T',
	'T':'A',
	'C':'G',
	'G':'C',
	'N':'N',
	'R':'Y',
	'Y':'R',
	'W':'W',
	'S':'S',
	'M':'K',
	'K':'M',
	'B':'V',
	'V':'B',
	'H':'D',
	'D':'H'}
	nucleotides = list(sequence)
	complement = [complements_dict[nucleotide] for nucleotide in nucleotides]
	complement = ''.join(complement)
	reverse_complement = complement[::-1]
	return reverse_complement

def FindGuides(parsed_pam_df, snp_df, chr_dict):
	# Takes in a parsed_pam_df and a snp_df.
	# Find guides for each Cas in the genomic context of each SNP.
	# Outputs a two-element tuple with a sgrna_df (indexed by sgRNA) and a summary_df (indexed by SNP).

	pam_location_list = []
	summary_list = []
	# Create a list of column headers for parsed_pam_df, FindGuides input values, and sgrna_df.
	snp_df_columns_list = list(snp_df.columns.values)
	pam_df_columns_list = list(parsed_pam_df.columns.values)
	input_values_labels = snp_df_columns_list[:]
	input_values_labels.extend(pam_df_columns_list)
	input_values_labels.append('Strand')
	sgrna_df_columns_list = input_values_labels[:]
	sgrna_df_columns_list.extend(['PAM_Window','sgRNA','PAM_Location','Found_PAM','Editing_Window','Num_Alt_Alleles'])

	for i in range(0,len(snp_df)):
		print('Finding guides for SNP # ' + str(i))
		r = list(snp_df.iloc[i,:])
		# pam_count_list stores the pam counts for each row.
		pam_count_list = r[:]
		# The pam counts for each Cas are appended to temp_pam_count_list.
		temp_pam_count_list = []
		# Create a temporary dictionary to map position in sequence to position in ref genome.
		sequence = r[snp_df_columns_list.index('Genomic_Context')]
		temp_genomic_location_dict = {}
		for base in range(0,len(sequence)):
			temp_genomic_location_dict[base] = r[snp_df_columns_list.index('Genomic_Context_Start')] + base
		# Go through each Cas in parsed_pam_df.
		for index in range(0,len(parsed_pam_df)):
			row = list(parsed_pam_df.iloc[index,:])
			# Initialize PAM count to 0.
			pam_count = 0
			# Create a list of input values for FindGuides in sense strand.
			input_values = r[:]
			input_values.extend(row)
			# Search for guides in sense strand.
			input_values_sense = input_values[:]
			input_values_sense.append('sense')
			find_guides_output_sense = FindGuidesInSequence(input_values_sense, input_values_labels, chr_dict, temp_genomic_location_dict)
			for guide in range(0,len(find_guides_output_sense)):
				# If a guide is found, increment the PAM count.
				if find_guides_output_sense[guide][sgrna_df_columns_list.index('PAM_Location')] != 'NaN':
					pam_count += 1
				pam_location_list.append(find_guides_output_sense[guide])
			# Create a list of input values for FindGuides in antisense strand.
			input_values_antisense = input_values[:]
			input_values_antisense.append('antisense')
			input_values_antisense[input_values_labels.index('Genomic_Context')] = ReverseComplement(r[snp_df_columns_list.index('Genomic_Context')])
			input_values_antisense[input_values_labels.index('ReferenceAllele')] = ReverseComplement(r[snp_df_columns_list.index('ReferenceAllele')])
			input_values_antisense[input_values_labels.index('AlternateAllele')] = ReverseComplement(r[snp_df_columns_list.index('AlternateAllele')])
			# Search for guides in antisense strand.
			find_guides_output_antisense = FindGuidesInSequence(input_values_antisense, input_values_labels, chr_dict, temp_genomic_location_dict)
			for guide in range(0,len(find_guides_output_antisense)):
				# Reset ref allele, alt allele, and genomic context.
				find_guides_output_antisense[guide][sgrna_df_columns_list.index('ReferenceAllele')] = r[snp_df_columns_list.index('ReferenceAllele')]
				find_guides_output_antisense[guide][sgrna_df_columns_list.index('AlternateAllele')] = r[snp_df_columns_list.index('AlternateAllele')]
				find_guides_output_antisense[guide][sgrna_df_columns_list.index('Genomic_Context')] = r[snp_df_columns_list.index('Genomic_Context')]
				# If a guide is found, reset PAM_Location (still indicating the 5' end of PAM) and increment pam_count.
				if find_guides_output_antisense[guide][sgrna_df_columns_list.index('PAM_Location')] != 'NaN':
					find_guides_output_antisense[guide][sgrna_df_columns_list.index('PAM_Location')] = len(sequence)-1-find_guides_output_antisense[guide][sgrna_df_columns_list.index('PAM_Location')]
					pam_count += 1
				pam_location_list.append(find_guides_output_antisense[guide])
			temp_pam_count_list.append(pam_count)
		# Get total number of PAM sites.
		total_pam_count = sum(temp_pam_count_list)
		temp_pam_count_list.append(total_pam_count)
		pam_count_list.extend(temp_pam_count_list)
		summary_list.append(pam_count_list)
	# Create summary dataframe indexed by SNP.
	summary_df_columns_list = snp_df_columns_list[:]
	summary_df_columns_list.extend(parsed_pam_df['Cas'])
	summary_df_columns_list.extend(['PAM_Count'])	
	if summary_list:
		summary_df = pd.DataFrame(data=summary_list, columns=summary_df_columns_list)
	if not summary_list:
		print('Error: summary_list is empty')
		return
	# Create sgRNA dataframe indexed by sgRNA.
	if pam_location_list:
		sgrna_df = pd.DataFrame(data=pam_location_list, columns=sgrna_df_columns_list)
	if not pam_location_list:
		print('Error: pam_location_list is empty')
		return
	return (sgrna_df, summary_df)

def FindGuidesInSequence(input_values, input_values_labels, chr_dict, temp_genomic_location_dict):
	# Takes in information about SNP, sequence, and PAM.
	# Finds guides using CountPAMs.
	# Outputs a list of guides in given sequence.
	
	nucleotide_codes_dict = {'A':{'A'},
	'C':{'C'},
	'G':{'G'},
	'T':{'T'},
	'R':{'A','G'},
	'Y':{'C','T'},
	'W':{'A','T'},
	'S':{'G','C'},
	'M':{'A','C'},
	'K':{'G','T'},
	'B':{'G','C','T','Y','S','K'},
	'H':{'A','C','T','Y','W','M'},
	'D':{'A','G','T','R','W','K'},
	'V':{'A','G','C','R','S','M'},
	'N':{'A','G','C','T','R','Y','W','S','M','K','B','H','D','V'},}

	sgrna_list = []
	allowed_ref_allele = input_values[input_values_labels.index('Allowed_Ref_Allele')]
	allowed_alt_allele = input_values[input_values_labels.index('Allowed_Alt_Allele')]
	allowed_ref_allele_set = nucleotide_codes_dict[allowed_ref_allele]
	allowed_alt_allele_set = nucleotide_codes_dict[allowed_alt_allele]
	# Inititalize pam_window, sgRNA, pam_location_in_sequence, found_pam, and num_alt_alleles to null values.
	pam_window = 'NaN'
	sgrna = 'No sgRNA found'
	pam_location_in_sequence = 'NaN'
	found_pam = 'NaN'
	num_alt_alleles = 'NaN'
	editing_window = 'NaN'
	# Check whether nucleotide is allowed.
	if input_values[input_values_labels.index('ReferenceAllele')] not in allowed_ref_allele_set:
		input_values.extend([pam_window,sgrna,pam_location_in_sequence,found_pam,editing_window,num_alt_alleles])
		sgrna_list.append(input_values)
		return sgrna_list
	if input_values[input_values_labels.index('AlternateAllele')] not in allowed_alt_allele_set:
		input_values.extend([pam_window,sgrna,pam_location_in_sequence,found_pam,editing_window,num_alt_alleles])
		sgrna_list.append(input_values)
		return sgrna_list
	# Get input values.
	sequence = input_values[input_values_labels.index('Genomic_Context')]
	pam_window_start = input_values[input_values_labels.index('PAM_Window_Start')]
	pam_window_end = input_values[input_values_labels.index('PAM_Window_End')]
	# Calculate snp_position_in_sequence (should be 50).
	snp_position_in_sequence = input_values[input_values_labels.index('SNP_Position')] - input_values[input_values_labels.index('Genomic_Context_Start')]
	# Go through each PAM in the list of input PAMs.
	pam_list = input_values[input_values_labels.index('Split_PAM')]
	pamside = input_values[input_values_labels.index('PAMside')]
	strand = input_values[input_values_labels.index('Strand')]
	for pam in pam_list:
		# Slice out pam window, taking into account the length of the PAM.
		if pamside == 3:
			pam_window = sequence[(snp_position_in_sequence+pam_window_start):(snp_position_in_sequence+pam_window_end+len(pam))]
		if pamside == 5:
			pam_window = sequence[(snp_position_in_sequence+pam_window_start-(len(pam)-1)):(snp_position_in_sequence+pam_window_end+1)]
		# Check whether given genomic context is long enough for pam window.
		if snp_position_in_sequence+pam_window_start < 0 or snp_position_in_sequence+pam_window_end+len(pam) > len(sequence)-1:
			sgrna = 'Error: genomic context is not long enough for the pam window'
			input_values.extend([pam_window,sgrna,pam_location_in_sequence,found_pam,editing_window,num_alt_alleles])
			sgrna_list.append(input_values)
			return sgrna_list
		# Search for pam in strand.
		# CountPAMs returns pam_location_in_pam_window.
		pams_in_pam_window_list = CountPAMs(pam_window,pam)
		if not pams_in_pam_window_list:
			sgrna = 'No sgRNA found'
			pam_location_in_sequence = 'NaN'
			found_pam = 'NaN'
			num_alt_alleles = 'NaN'
			editing_window = 'NaN'
			temp_input_values = input_values[:]
			temp_input_values.extend([pam_window,sgrna,pam_location_in_sequence,found_pam,editing_window,num_alt_alleles])
			sgrna_list.append(temp_input_values)
		if pams_in_pam_window_list:
			for pam_location in pams_in_pam_window_list:
				temp_input_values = input_values[:]
				if pamside == 3:
					pam_location_in_sequence = pam_location + snp_position_in_sequence + pam_window_start
					# Check whether the alt allele is unique in the editing window (important for base editors).
					editing_window = sequence[(pam_location_in_sequence-pam_window_end):(pam_location_in_sequence-pam_window_start)+1]
					if allowed_alt_allele != 'N':
						num_alt_alleles = editing_window.count(allowed_alt_allele)
				if pamside == 5:
					# Subtract len(pam) because the pam_window is extended in the 5' direction for 5' PAMs.
					pam_location_in_sequence = pam_location + snp_position_in_sequence + pam_window_start - (len(pam)-1)
					# Check whether the alt allele is unique in the editing window.
					editing_window = sequence[(pam_location_in_sequence-pam_window_end+len(pam)):(pam_location_in_sequence-pam_window_start+len(pam)+1)]
					if allowed_alt_allele != 'N':
						num_alt_alleles = editing_window.count(allowed_alt_allele)
				# Extract guides.
				chrom_num = input_values[input_values_labels.index('Chromosome')]
				spacer_len = input_values[input_values_labels.index('Spacer_Len')]
				if pamside == 3:
					found_pam = sequence[pam_location_in_sequence:(pam_location_in_sequence+len(pam))]
				# If there is not enough genomic context for the spacer, go back and pull out the sgRNA sequence from the reference genome.	
					if pam_location_in_sequence-spacer_len < 0 or pam_location_in_sequence-spacer_len > len(sequence)-1:
						chr_sequence = chr_dict[chrom_num]
						sgrna = chr_sequence[(pam_location_in_sequence-spacer_len):pam_location_in_sequence]
					else:
						sgrna = sequence[(pam_location_in_sequence-spacer_len):pam_location_in_sequence]
				if pamside == 5:
					found_pam = sequence[pam_location_in_sequence:pam_location_in_sequence+len(pam)]
					if pam_location_in_sequence+len(pam) < 0 or pam_location_in_sequence+len(pam)+spacer_len > len(sequence)-1:
						chr_sequence = chr_dict[chrom_num]
						sgrna = chr_sequence[(pam_location_in_sequence+len(pam)):(pam_location_in_sequence+len(pam)+spacer_len)]
					else:
						sgrna = sequence[(pam_location_in_sequence+len(pam)):(pam_location_in_sequence+len(pam)+spacer_len)]
				temp_input_values.extend([pam_window,sgrna,pam_location_in_sequence,found_pam,editing_window,num_alt_alleles])
				sgrna_list.append(temp_input_values)
	return sgrna_list

def CountPAMs(sequence,pam):
	# Takes in sequence (5' to 3') and a PAM to search for.
	# Outputs a list of PAM locations, where the location indicates the 5' end of the PAM.
	# CountPAMs is agnostic to strand.
	# To find PAMs in the antisense strand, CountPAMs takes in the reverse complement (5' to 3') and outputs the 5' location of the PAM in the antisense strand.

	nucleotide_codes_dict = {'A':{'A'},
	'C':{'C'},
	'G':{'G'},
	'T':{'T'},
	'R':{'A','G'},
	'Y':{'C','T'},
	'W':{'A','T'},
	'S':{'G','C'},
	'M':{'A','C'},
	'K':{'G','T'},
	'B':{'G','C','T'},
	'H':{'A','C','T'},
	'D':{'A','G','T'},
	'V':{'A','G','C'},
	'N':{'A','G','C','T'},}

	pam_location_list = []
	for base in range(0, len(sequence)+1-len(pam)):
		for i in range (0, len(pam)):
			nucleotide = sequence[base+i]
			pam_allowed_nucleotides = nucleotide_codes_dict[pam[i]]
			pam_found = True
			if nucleotide not in pam_allowed_nucleotides:
				pam_found = False
				break
		if pam_found == True:
			pam_location_list.append(base)
	return pam_location_list

def main():
	# Read in and parse variant file from ClinVar.
	args = GetParser().parse_args()
	variant_file = args.variant_file
	variant_df = pd.read_table(variant_file, dtype = {'#AlleleID':str,'Name':str,'ClinicalSignificance':str,'Assembly':str,'Chromosome':str,'Type':str,'Start':int,'ReferenceAllele':str,'AlternateAllele':str}, usecols = ['#AlleleID','Name','ClinicalSignificance','Assembly','Chromosome','Type','Start','ReferenceAllele','AlternateAllele'])
	# Create output folder.
	output_name = args.output_name
	output_folder = output_name + '_' + str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	# Parse variant_df.
	parsed_variant_df = ParseVariantFile(variant_df)
	parsed_variant_df.to_csv(output_folder + '/parsed_variant_df.txt',sep='\t',index=False)
	# Read in chrom files.
	chr_dict = ReadChromFiles(1,24)
	# locate_snp_output is a two-element tuple.
	locate_snp_output = LocateSNP(chr_dict, parsed_variant_df)
	snp_error_df = locate_snp_output[0]
	if not snp_error_df.empty:
		snp_error_df.to_csv(output_folder + '/snp_error_df.txt',sep='\t',index=False)
	snp_df = locate_snp_output[1]
	snp_df.to_csv(output_folder + '/snp_df.txt',sep='\t',index=False)
	# Read in pam file and parse list of PAMs into list.
	pam_file = args.pam_file
	pam_df = pd.read_table(pam_file)
	parsed_pam_df = ParsePAMFile(pam_df)
	# Count pams in each snp and produce list of sgRNAs.
	# find_guides_output is a two-element tuple.
	find_guides_output = FindGuides(parsed_pam_df,snp_df,chr_dict)
	sgrna_df = find_guides_output[0]
	sgrna_df.to_csv(output_folder + '/sgrna_df.txt', sep='\t',index=False)
	summary_df = find_guides_output[1]
	summary_df.to_csv(output_folder + '/summary_df.txt', sep='\t',index=False)
	# Create README.
	with open(output_folder+'/README.txt','w') as o:
		w = csv.writer(o,delimiter='\t')
		w.writerow((['Code Version: 1.12']))
		w.writerow((['Variant file:'+variant_file]))
		w.writerow((['PAM file:'+pam_file]))
		w.writerow((['Output folder:'+output_folder]))

if __name__ == '__main__':
	main()