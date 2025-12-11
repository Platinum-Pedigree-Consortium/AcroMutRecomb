import pysam
import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm


def get_bam(sample, parent, role, region, bams, sequence_platform):
	"""
	Retrieves the path to the sample's bam file.
	:param sample: De novo child (unique ID).
	:param parent: Reference parent (father or mother).
	:param role: Person to find bam for (child or parent).
	:param region: Region name.
	:param bams: Manifest dataframe of bam files.
	:param sequence platform: illumina, ont, or hifi
	:returns: PySam AlignmentFile object of the sample's read data.
	"""
	if role in ['father', 'mother']:
		person = 'parent'
	else:
		person = 'child'

	bam_path =  bams.loc[(bams['sample'] == sample) & (bams['parent_type'] == parent) & (bams['region'] == region)][f'{person}_{sequence_platform}'].values[0]
	if not bam_path:
		return None
	bam = pysam.AlignmentFile(bam_path, "rc")
	return bam

def get_mean_qual(qualities):
	"""
	Computes the mean base quality for reads.
	:param qualities: List of integer base qualities from read data.
	:return: Number of quality values in list and their mean.
	"""
	qual_count = len(qualities)

	if qual_count != 0:
		return qual_count, np.mean(qualities)

	return qual_count, "NA"

def check_snv_read(read, position, ref, alt, ref_quality, alt_quality, other_quality, map_quality):
	"""
	Evaluate a read to record its quality data and the base at the query position.
	:param read: AlignedSegment object of read aligned to query position.
	:param position: Integer query position of DNM in reference genome.
	:param ref: Reference allele, A/C/G/T.
	:param alt: de novo alternate allele, A/C/G/T.
	:param ref_quality: List of base quality values for reads with reference allele at query position.
	:param alt_quality: List of base quality values for reads with alternate allele at query position.
	:param other_quality: List of base quality values for reads with other allele at query position.
	:param map_quality: List of map quality values for reads at query position.
	:returns: Name of the read, base and its quality at query position, median base quality across the read, length of the read, and its map quality.
	"""

	read_length = read.query_length
	mapq = read.mapping_quality
	map_quality.append(mapq)
	try:
		ref_idx = read.get_reference_positions(full_length=True).index(position - 1)
	except:
		return read.query_name, "deletion", "NA", "NA", read_length, mapq

	base = read.query_sequence[ref_idx]
	try:
		quality = read.query_qualities[ref_idx]
	except TypeError:
		quality = 0
	read_pos = ref_idx/read_length

	if base == ref:
		correct_quality = ref_quality
	elif base == alt:
		correct_quality = alt_quality
	else:
		correct_quality = other_quality

	correct_quality.append(quality)

	return read.query_name, base, quality, read_pos, read_length, mapq

def get_snv_data(bam, site_id, chrom, pos, ref, alt, reads_out, samples_out):
	"""
	Write individual and summary stats to files for all the reads aligned to a query position in one sample.
	:param bam: AlignmentFile object of all reads from a sample's BAM.
	:site_id: Unique ID for DNM site, with sample name, chromosome, position, ref, and alt alleles.
	:chrom: Chromosome name for DNM site.
	:pos: Integer reference position of DNM site.
	:param ref: Reference allele, A/C/G/T.
	:param alt: de novo alternate allele, A/C/G/T.
	:param reads_out: Opened file of read data.
	:param samples_out: Opened file of summary data for each sample.
	"""
	sample_name, bam = bam
	ref_qual = []
	alt_qual = []
	other_qual = []
	map_qual = []

	for aligned_segment in bam.fetch(chrom, pos - 1, pos):
		if not aligned_segment.query_sequence:
			continue

		line = "\t".join((
			site_id,
			sample_name,
			*map(str, check_snv_read(aligned_segment, pos, ref, alt, ref_qual, alt_qual, other_qual, map_qual)),
		))
		reads_out.write(line)
		reads_out.write("\n")


	ref_count, avg_ref_q = get_mean_qual(ref_qual)
	alt_count, avg_alt_q = get_mean_qual(alt_qual)
	other_count, avg_other_q = get_mean_qual(other_qual)
	read_count, avg_mapq = get_mean_qual(map_qual)

	sample_line = "\t".join(map(str, (site_id, sample_name,
							ref_count, avg_ref_q,
							alt_count, avg_alt_q,
							other_count, avg_other_q,
							avg_mapq)))
	samples_out.write(sample_line)
	samples_out.write("\n")


def main(variants, bams, platform, denovo_sample, parent, region):
	variant_df = pd.read_csv(variants, sep='\t', header = None)

	bam_df = pd.read_csv(bams, sep='\t', dtype=str, keep_default_na=False)
	bam_df = bam_df.astype(str)

	sample_dict = {parent:  bam_df.loc[(bam_df['sample'] == denovo_sample) & (bam_df['parent_type'] == parent)]['parent_id'].values[0],
				   'child': denovo_sample}


	reads_out = open(f'read_data/{denovo_sample}_from_{parent}.{region}.{platform}.read_data.tsv', 'w')
	samples_out = open(f'read_data/{denovo_sample}_from_{parent}.{region}.{platform}.read_summary.tsv', 'w')

	# write headers to output files

	print("\t".join(("id", "sample", "read_id", "allele", "base_quality", "read_position",
				"read_length", "mapq")), file = reads_out)
	print("\t".join(("id", "sample", "ref_count", "mean_ref_qual", "alt_count", "mean_alt_qual",
				"other_count", "mean_other_qual", "avg_mapq")), file = samples_out)

	for role, sample_id in sample_dict.items():
		print(f"currently working on sample {sample_id}, now processing:")
		sample_bam = (sample_id, get_bam(denovo_sample, parent, role, region, bam_df, platform))
		if not sample_bam[1]:
			continue

		for idx, row in tqdm(variant_df.iterrows(), total=len(variant_df)):
			if row[0] == 'chr':
				continue
			get_snv_data(sample_bam, row[2], row[0], int(row[1]), row[3], row[4], reads_out, samples_out)

	reads_out.close()
	samples_out.close()

if __name__ == "__main__":
	# to run from the command line, use commented out lines
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	# ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
	# ap.add_argument("-c", "--child", required=True, type=str, help="sample id")
	# ap.add_argument("-f", "--parent", required=True, type=str, help="parent id")
	# ap.add_argument("-r", "--region", required=True, type=str, help="region")
	# ap.add_argument("-p", "--platform", required=True, help="sequencing platform")
	#
	#
	# args = ap.parse_args()
	#
	# main(args.input, args.bams, args.platform, args.child, args.parent, args.region)
	main(snakemake.input.snv_calls, snakemake.params.bam_manifest, snakemake.params.platform, snakemake.params.child, snakemake.params.parent, snakemake.params.region)
