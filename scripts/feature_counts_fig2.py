import sys, HTSeq, argparse
from collections import defaultdict

def main(args):
 
	premature = defaultdict(int)
	mature = defaultdict(int)
	introns = set()

	features = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	for line in open(args.feature_ranges):
		fields = line.rstrip().split('\t')
		iv = HTSeq.GenomicInterval(fields[0], int(fields[1])-1, int(fields[2]), fields[4])
		features[iv] += fields[3]
		introns.add(fields[3]) 

	for alignment in HTSeq.BAM_Reader(args.input_file):
		
		feat = set()
		junctions = set()

		for cigop in alignment.cigar:
			if cigop.type == 'M':
				for iv, val in features[cigop.ref_iv].steps():
					feat |= val
			elif cigop.type == 'N':
				for iv, val in features[cigop.ref_iv].steps():
					junctions |= val


		candidate_genes = set()
		for i in feat:
			candidate_genes.add(i.split(';')[0])
		for i in junctions:
			candidate_genes.add(i.split(';')[0])

		if len(candidate_genes) == 1:

			for i in junctions:
				mature[i] += 1
			for i in feat:
				premature[i] += 1			



	print 'Intron\tMature\tPremature'
	for intron in sorted(introns):
		print '%s\t%d\t%d' % (intron, mature[intron], premature[intron])


def parseArguments():
	parser = argparse.ArgumentParser(prog="feature_counts_fig2.py", description='Counts the number of premature and mature reads in a .bam file.', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', '--input_file', required=True, help=' Input alignment file. (.bam)', metavar='', dest='input_file')
	required.add_argument('-f', '--features', required=True, help=' Feature ranges. (.bed)', metavar='', dest='feature_ranges')
	return parser.parse_args()

args = parseArguments()
main(args)
