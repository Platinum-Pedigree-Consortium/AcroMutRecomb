import argparse
import numpy as np
from typing import Tuple, List, Dict, Iterable, Optional
from collections import Counter
from collections import defaultdict


def get_var_info(var: List[str]):
	"""
	Subset relevant variant info for the output file.
	:param var: List of variant data from the input variant file.
	:return: List of relevant info, including variant identifiers and caller.
	"""
	var_info = var[0:5]
	return var_info


def import_reads(read_file: str):
	"""
	Import sample's read_data file to a dictionary of sample read data mapped to variants.
	:param read_file: Path to read_data file.
	:return: Dict of dicts, for each variant key, there should be a dict of read data mapped to each family member.
	"""
	read_dict = defaultdict(lambda: defaultdict(list))

	with open(read_file, 'r') as reads:
		header = reads.readline().rstrip().split('\t')
		basequal_idx = header.index('base_quality')
		mapq_idx = header.index('mapq')

		for line in reads:
			read = line.rstrip().split('\t')
			var_id, sample = read[0:2]

			read_data = (read[2], read[3], read[basequal_idx], read[mapq_idx])
			read_dict[var_id][sample].append(read_data)

	return read_dict

def mapq_filter(mapq: str, threshold: int):
	"""
	Filter mapq based on threshold.
	:param mapq: Mapq value for a given read.
	:param threshold: Mapq threshold as defined for sequencing platform.
	:return: True if above threshold, False if below threshold.
	"""
	try:
		 mapq = int(mapq)
	except ValueError:
		 return False

	if mapq < threshold:
		return False

	return True



def filter_parent_alleles(read_list: List[str]):
	"""
	Filter parental reads based on mapq, and then assign alleles based on base quality.
	:param read_list: List of read data, including read name, allele, base quality, and mapq.
	:return: List of alleles with base quality under 20, and alleles with base quality 20 or above.
	"""
	mapq_threshold = mapq_thresholds[platform]
	alleles = []

	for read in read_list:
		if not mapq_filter(read[3], mapq_threshold):
			continue

		alleles.append(read[1])

	return alleles



def filter_child_alleles(read_list: List[str]):
	"""
	Filter parental reads based on mapq, and then assign alleles based on base quality.
	:param read_list: List of read data, including read name, allele, base quality, and mapq.
	:return: List of alleles from blood, and alleles from cell line.
	"""
	mapq_threshold = mapq_thresholds[platform]
	allele_list = []

	for read in read_list:
		if not mapq_filter(read[3], mapq_threshold):
			continue
		allele_list.append(read[1])

	return allele_list


def check_allele_count_threshold(allele_list: list, quality: str, alt_allele: str):
	"""
	Determine if alternate allele may have been inherited from list of parental alleles.
	:param allele_list: Either high or low quality parental allele list.
	:param quality: 'high' if base quality 20 or greater, 'low' if not.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: True if alternate allele count is below threshold, False if it is equal or above threshold.
	"""
	threshold = inheritance_thresholds[platform][quality]
	alt_allele_count = Counter(allele_list)[alt_allele]

	if alt_allele_count > threshold:
		return False

	return True


def determine_possible_inheritance(alleles: list,  alt_allele: str):
	"""
	Count the high and low quality alternate alleles, and compare them to inheritance threshold.
	:param high_qual_alleles: List of alleles with base quality 20 or greater.
	:param low_qual_alleles: List of alleles with base quality under 20.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: Tuple of True/False if total allele count is greater than 0, True/False if high qual alts below threshold, True/False if low qual alts below threshold.
	"""

	depth = len([x for x in alleles if x != 'deletion'])

	if depth < 5:
		return (False, False)

	high_qual_pass = check_allele_count_threshold(alleles, 'high', alt_allele)

	return (True, high_qual_pass)


def calculate_depth_and_ab(allele_list: list, alt_allele: str):
	"""
	Calculate read depth and allele balance for blood or cell data.
	:param allele_list: List of alleles.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: Tuple of alt allele count, total read depth, and alt allele balance.
	"""
	depth = len(allele_list)
	if depth == 0:
		return (0, 0, np.nan)

	alt_count = Counter(allele_list)[alt_allele]
	ab = round(alt_count / depth, 5)

	return (alt_count, depth, ab)


def determine_false_positive(allele_list: list, allele_stats: Tuple):
	"""
	Determine if variant is true, inherited, false, or missing data.
	:param allele_list: List of alleles.
	:param allele_stats: Tuple with child alt allele count, read depth, and allele balance.
	:return: Tuple of True/False if depth is greater than 0, True/False if alt allele is represented.
	"""
	if len([x for x in allele_list if x != 'deletion']) < 5:
		return (False, False)

	alt_threshold = 0
	if platform == 'hifi':
		alt_threshold = 1

	if allele_stats[0] > alt_threshold:
		return (True, True)

	return (True, False)


def validate_variant(child_summary: Tuple, parent_summary: Tuple):
	"""
	Determine if variant is true, inherited, false, or missing data.
	:param chidl_summary: Tuple of True/False if child has blood data, True/False if child has alt allele in blood.
	:param parent_summary: Tuple of True/False if dad has data, True/False if alt in high qual reads, True/False if alt in low qual reads.
	:return: Final validation of variant.
	"""
	if parent_summary[0] == False:
		return 'missing_parent_data'

	if False in parent_summary[1:]:
		return 'inherited'

	if child_summary[0] == False:
		return 'missing_child_data'

	if child_summary[1] == True:
		return 'true_de_novo'

	return 'false_positive'


def get_parent_data(variant_reads, alt_allele):
	"""
	Count the number of alternate alleles on filtered reads and determine possible inheritance.
	:param variant_read_dict: List of reads containing variant from parent.
	:param alt_allele: Alternate allele from de novo mutation.
	:returns:  Tuple of True/False if parent has data, True/False if alt in high qual reads, True/False if alt in low qual reads.
	"""
	parent_alleles = filter_parent_alleles(variant_reads)
	return determine_possible_inheritance(parent_alleles, alt_allele)


def main(variant_file, read_file, parent, child, output_file):
	out_header = ['chr', 'pos', 'id', 'ref', 'alt', 'child_filtered_count', 'child_depth', 'child_allele_balance', 'filtered_status']

	variant_reads = import_reads(read_file)

	with open(variant_file, 'r') as variants, open(output_file, 'w') as outfile:
		print('\t'.join(out_header), file = outfile)
		for line in variants:
			var_data = line.rstrip().split('\t')
			id, ref, alt = var_data[2:5]

			variant_info = get_var_info(var_data)

			parent_inh_summary = get_parent_data(variant_reads[id][parent], alt)

			child_alleles = filter_child_alleles(variant_reads[id][child])
			child_stats = calculate_depth_and_ab(child_alleles, alt)
			child_dnm = determine_false_positive(child_alleles, child_stats)

			validation = validate_variant(child_dnm, parent_inh_summary)

			out_data = variant_info + [str(x) for x in child_stats] + [validation]

			print('\t'.join(out_data), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="Candidate de novo file")
	# ap.add_argument("-r", "--reads", required=True, help="File of read data")
	# ap.add_argument("-p", "--parent", required=True, help="parent ID", type=str)
	# ap.add_argument("-c", "--child", required=True, help="sample ID", type=str)
	# ap.add_argument("-o", "--output", required=True, help="Output file")
	# ap.add_argument("-p", "--platform", required=True, help="Output file")
	# args = ap.parse_args()
	#
	# platform = args.platform
	platform = snakemake.params.platform

	inheritance_thresholds = {
		'hifi': {'high': 1, 'low': 2},
		'ont': {'high': 2, 'low': 3},
		'illumina': {'high': 1, 'low': 2}
	}

	mapq_thresholds = {
		'hifi': 0,
		'ont': 0,
		'illumina': 0
	}

	# main(args.variants, args.reads,  args.parent,  args.child, args.output)
	main(snakemake.input.snv_calls,
		snakemake.input.read_data,
		snakemake.params.parent,
		snakemake.params.child,
		snakemake.output.platform_validation)
