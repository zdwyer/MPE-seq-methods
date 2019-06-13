import sys, HTSeq, argparse
from collections import defaultdict

def main(args):
 
	counts = {}
	problem = 0

	features = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	for line in open(args.feature_ranges):
		fields = line.rstrip().split('\t')
		iv = HTSeq.GenomicInterval(fields[0], int(fields[1])-1, int(fields[2]), fields[4])
		features[iv] += fields[3]
		counts[fields[3].split(';')[0]] = defaultdict(int)

	for alignment in HTSeq.BAM_Reader(args.input_file):
		

		feat = set()
		junctions = set()

		for i, cigop in enumerate(alignment.cigar):
			if cigop.type == 'M':
				for iv, val in features[cigop.ref_iv].steps():
					feat |= val
			elif cigop.type == 'N':
				if alignment.cigar[i-1].type =='M' and alignment.cigar[i-1].size > 3 and alignment.cigar[i+1].type =='M' and alignment.cigar[i+1].size > 3:
					for iv, val in features[cigop.ref_iv].steps():
						junctions |= val

		candidate_genes = set()
		for i in feat:
			candidate_genes.add(i.split(';')[0])
		for i in junctions:
			candidate_genes.add(i.split(';')[0])

		if len(candidate_genes) == 1:
			if len(junctions) == 1:
				counts[list(junctions)[0].split(';')[0]]['Mature'] += 1
			else:
				if len(feat) > 0:
					gene = list(feat)[0].split(';')[0]
					pieces = set()
					for i in feat:
						pieces.add(i.split(';')[1])
					if 'I1' in pieces and ('E1' in pieces or 'E2' in pieces):
						counts[gene]['Premature'] += 1
					elif 'E1' in pieces or 'E2' in pieces:
						counts[gene]['Exonic'] += 1

	print 'Gene\tExonic\tMature\tPremature'
	for gene in sorted(counts):
		print '%s\t%d\t%d\t%d' % (gene, counts[gene]['Exonic'], counts[gene]['Mature'], counts[gene]['Premature'])

def parseArguments():
	parser = argparse.ArgumentParser(prog="", description='', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', '--input_file', required=True, help=' Input alignment file. (.bam)', metavar='', dest='input_file')
	required.add_argument('-f', '--features', required=True, help=' Feature ranges. (.bed)', metavar='', dest='feature_ranges')
	return parser.parse_args()

args = parseArguments()
main(args)