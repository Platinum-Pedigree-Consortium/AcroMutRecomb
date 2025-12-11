import argparse

def get_platform(file_name):
	"""
	Infer platform name from file name.
	:param file_name: Path to platform validation file.
	:return platform: 'hifi' 'ont' or 'illumina'
	"""
	if 'hifi' in file_name:
		return 'hifi'
	elif 'ont' in file_name:
		return 'ont'
		

def add_platform(platform, var_file):
	"""
	Add platform to the start of every line in corresponding validation file.
	:param platform: 'hifi' 'ont' or 'illumina'
	:param var_fil: Opened validation file.
	:return: Lines from file wtiht platform at the start.
	"""
	next(var_file)
	for line in var_file:
		yield platform, line.rstrip().split('\t')

def cross_platform_validated(pf_data):
	"""
	Determines final validation based on each of the three sequencing platforms.
	:param pf_data: Dictionary of platform validations mapped to platform names.
	:returns: Final validation (eg 'true_de_novo', 'inherited').
	"""
	hifi_val = pf_data['hifi']
	ont_val = pf_data['ont']

	if hifi_val == 'true_de_novo' and ont_val == 'true_de_novo':
		return 'true_de_novo'

	if hifi_val == 'true_de_novo' and ont_val == 'missing_parent_data':
		return 'hifi_true_de_novo_ont_no_parent_data'

	if hifi_val == 'true_de_novo' and ont_val == 'false_positive':
		return 'hifi_true_de_novo_ont_false_positive'

	if 'inherited' in [hifi_val, ont_val]:
		return 'inherited'

	if hifi_val == 'false_positive':
		return 'false_positive'

	if 'missing' in hifi_val:
		if 'missing_data' in [ont_val]:
			return 'cannot_validate'
		else:
			return 'false_positive'

	return f'missing_case_hifi_{hifi_val}_ont_{ont_val}'

def validate_variant(variant_data):
	"""
	Makes dictionary of platform validations and runs final validation script.
	:param variant_data: List of lists corresponding to each validation file.
	:return: Final validation (eg 'true_de_novo', 'inherited').
	"""
	platform_validation = {x[0]: x[1][-1] for x in variant_data}
	final_validation = cross_platform_validated(platform_validation)

	return final_validation

def main(validation_files, output_file):
	variant_files = {get_platform(path): path for path in validation_files}
	files_with_platform = [add_platform(pf, open(path, 'r')) for pf, path in variant_files.items()]

	out_columns = ['chr', 'pos', 'id', 'ref', 'alt'] + \
		[f'{x}_alt_allele_count' for x in  variant_files.keys()] + \
		[f'{x}_dp' for x in  variant_files.keys()] + \
		[f'{x}_ab' for x in  variant_files.keys()] + \
		['final_validation']

	with open(output_file, 'w') as outfile:
		print('\t'.join(out_columns), file=outfile)

		for variants in zip(*files_with_platform):
			variant = variants[0][1][:5]
			for platform, line in variants:
				if line[:5] != variant:
					raise ValueError("your files are not in sorted order!")
			val = validate_variant(variants)
			# if 'true_de_novo' in val:
			variant_out = variant + [x[1][-4] for x in variants] + [x[1][-3] for x in variants] + [x[1][-2] for x in variants] + [val]
			print('\t'.join(variant_out), file=outfile)


	for file in files_with_platform:
		try:
			next(file)
		except StopIteration:
			pass
		else:
			raise ValueError("one of your files is longer than the other!")

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-f", "--files", required=True, help="validation files", nargs='+')
	# ap.add_argument("-o", "--output", required=True, help="output")
	#
	# args = ap.parse_args()

	# main(args.files, args.output)
	main(snakemake.input.platform_validations, snakemake.output.validated_snvs)
