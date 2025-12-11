import argparse
import pysam
from pysam import AlignmentFile
from typing import Tuple, List, Dict

def import_contig_names(contig_file):
	contig_dict = {}
	with open(contig_file, 'r') as infile:
		for line in infile:
			sample, name = line.rstrip().split('\t')
			if sample not in contig_dict:
				contig_dict[sample] = set()

			contig_dict[sample].add(name)

	return contig_dict

def get_sample_information(idx, variant_info, contig_dict):
	sample_info = variant_info[idx].split('_')
	sample = sample_info[0]
	options = ['_'.join(sample_info[i:]) for i in range(len(sample_info))]
	chrom = next(contig for contig in options if contig in contig_dict[sample])

	pos = int(variant_info[idx+1])

	return sample, (chrom, pos)

def get_fasta_dict(manifest_file):
	fasta_dict = {}
	with open(manifest_file, 'r') as manifest:
		header = manifest.readline().rstrip().split('\t')
		fasta_idx = header.index('fasta')
		sm_idx = header.index('sample')

		for line in manifest:
			sample_info = line.rstrip().split('\t')
			fasta = pysam.FastaFile(sample_info[fasta_idx])
			fasta_dict[sample_info[sm_idx]] = fasta

	return fasta_dict

def reverse_complement(sequence):
	complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	complement = ''

	for base in sequence[::-1]:
		complement += complement_dict[base]

	return complement

class Variant:
	id: str
	parent: str
	parent_location: Tuple[str, int]
	child: str
	child_location: Tuple[str, int]
	strand: str
	ref: str
	alt: str
	parent_seq: str
	child_seq: str

	def __init__(self, data_line: str, fasta_dict: Dict[str, AlignmentFile], contig_dict: Dict[str, List[str]]):
		var_info = data_line.rstrip().split('\t')
		self.parent, self.parent_location = get_sample_information(0, var_info, contig_dict)
		self.child, self.child_location = get_sample_information(4, var_info, contig_dict)
		self.strand = var_info[3]

		self.id =  self.child + '_' + '_'.join(str(x) for x in self.child_location)

		self.parent_seq = self._get_allele(self.parent, self.parent_location, self.strand, fasta_dict)
		self.child_seq = self._get_allele(self.child, self.child_location, '+', fasta_dict)

		self.ref = self.parent_seq[4]
		self.alt = self.child_seq[4]


	def __repr__(self):
		return self.id

	def _get_allele(self, sample, location, strand, fasta_dict):
		chrom, pos = location
		fasta = fasta_dict[sample]

		if strand == '+':
			allele = fasta.fetch(chrom, pos - 5, pos + 4)
		else:
			allele = fasta.fetch(chrom, pos - 6, pos + 3)
			allele = reverse_complement(allele)

		return allele

	def print(self):
		out_data = [*self.child_location, self.id, self.ref, self.alt]
		out_data = [str(x) for x in out_data]
		return '\t'.join(out_data)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	ap.add_argument("-m", "--manifest", required=True, help="manifest file")
	ap.add_argument("-c", "--contig_names", required=True, help="contig names")
	args = ap.parse_args()

	contig_names = import_contig_names(args.contig_names)
	opened_fastas = get_fasta_dict(args.manifest)
	snv_list = []

	with open(args.input, 'r') as infile:
		header = infile.readline().rstrip().split('\t')
		for line in infile:
			snv = Variant(line, opened_fastas, contig_names)
			snv_list.append(snv)

			if snv.parent_seq[:4] != snv.child_seq[:4] or snv.parent_seq[5:] != snv.child_seq[5:]:
				continue

			print(snv.print())

