from pybedtools import BedTool
from pysam import FastaFile
from typing import Tuple, List, Dict, Iterable, Optional
from tqdm import tqdm
import argparse


def get_complements(complement_file: str):
	"""
	Read in triplet completment lookup table and save as dictionary.
	:param complement_file: Path to complement lookup table file.
	:return" Dictionary of triplets mapped to their complement.
	"""
	complements = {}
	with open(complement_file, 'r') as f:
		for line in f:
			key, val = line.rstrip().split('\t')
			complements[key] = val
			complements[val] = val

	return complements


def read_variant_file(variant_file: str, out_header = List[str]):
	"""
	Read in de novo variants, save them to a dictionary with bed intervals.
	:param variant_file: Path to file of de novo mutations.
	:return: Header of variant file as list of strings.
	:return: Dictionary of variant info mapped to bed interval.
	"""
	variants = {}
	with open(variant_file, 'r') as var_file:
		header = var_file.readline().rstrip().split('\t')
		for line in var_file:
			var_data = line.rstrip().split()
			pos = int(var_data[header.index('pos')])
			bed_info = (var_data[header.index('chr')], pos - 1, pos)
			if var_data[header.index('final_validation')] == 'true_de_novo':
				relevant_data = var_data[:5] + var_data[9:11]
				variants[bed_info] = {k:v for k,v in zip(out_header, relevant_data)}

	return variants


def transition_transversion(ref: str, alt: str):
	"""
	Assign variant as transition or transversion mutation.
	:param ref: Reference allele.
	:param alt: Alternate allele.
	:return: Transition or transversion.
	"""
	bases = set([ref, alt])
	if bases == set(['A', 'G']) or bases == set(['C', 'T']):
		return 'transition'
	else:
		return 'transversion'

def check_triplet(variant_info: List[str], ref_fasta: FastaFile, complement_lookup: Dict[str, str]):
	"""
	Create triplet name for variant, with 3' and 5' base converted into correct complement, and describes mutation type.
	:param variant_info: List of variant data from input file.
	:param ref_fasta: FastaFile object of reference genome sequence.
	:complement_lookup: Dictionary of complement sequence mapped to correct complement.
	:return: Dictionary of triplet sequence, whether the mutation is CpG>CpG, and substitution_type.
	"""
	ref = variant_info['ref']
	alt = variant_info['alt']
	pos1 = int(variant_info['pos']) - 2
	pos2 = int(variant_info['pos']) + 1
	seq = ref_fasta.fetch(variant_info['chr'], start=pos1, end=pos2).upper()
	triplet = "%s-%s"%(seq, alt)
	triplet = complement_lookup[triplet]

	if triplet[-4:-2] == 'CG':
		CpG = 'True'
	else:
		CpG = 'False'

	trans = transition_transversion(ref, alt)

	return {'tiplet': triplet, 'CpG': CpG, 'substitution_type': trans}

def main(variant_file, fasta_file, complement_lookup, output_file, sample, parent, region):

	out_header = ['sample', 'parent', 'region', 'chr', 'pos', 'id', 'ref', 'alt', 'hifi_ab', 'ont_ab', 'triplet', 'CpG', 'substitution_type']

	variant_dict = read_variant_file(variant_file, out_header[3:10])
	reference = FastaFile(fasta_file)
	complement_dict = get_complements(complement_lookup)

	with open(output_file, 'w') as outfile:
		print('\t'.join(out_header), file = outfile)

		for interval, var_info in tqdm(variant_dict.items(), total = len(variant_dict.items())):
			var_info.update(check_triplet(var_info, reference, complement_dict))
			out_info = [sample, parent, region] + list(var_info.values())

			print('\t'.join(out_info), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="Variant file")
	# ap.add_argument("-f", "--fasta", required=True, type=str, help="Reference genome fasta")
	# ap.add_argument("-l", "--lookup", required=True, type=str, help="triplet complement lookup table")
	# ap.add_argument("-o", "--output", required=True, help="output file")
	# ap.add_argument("-s", "--sample", required=True, help="sample")
	# ap.add_argument("-p", "--parent", required=True, help="parent")
	# ap.add_argument("-r", "--region", required=True, help="region")
	#
	#
	# args = ap.parse_args()
	#
	# main(args.input, args.fasta, args.lookup, args.output, args.sample, args.parent, args.region)
	main(snakemake.input.dnm_file, snakemake.input.fasta, snakemake.params.lookup_table, snakemake.output.annotated_snvs, snakemake.params.sample, snakemake.params.parent, snakemake.params.region)
