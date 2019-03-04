import argparse, gzip
from itertools import izip_longest

def main(args):

	occured = set() # Keeps track of read 1/read 2 sequence combinations that have already occured
	unique = [] # Keeps a single instance of each read 1/read 2 sequence combination

	with gzip.open(args.read_1) as i1, gzip.open(args.read_2) as i2:
		for line1, line2 in izip_longest(i1, i2):
			
			#Read one entry from read 1
			info1 = line1.rstrip().split(' ')[0]
			seq1 = next(i1).rstrip()
			extra1 = next(i1).rstrip()
			quality1 = next(i1).rstrip()

			#Read one entry from read 2
			info2 = line2.rstrip().split(' ')[0]
			seq2 = next(i2).rstrip()
			extra2 = next(i2).rstrip()
			quality2 = next(i2).rstrip()

			# If this read 1/read 2 sequence combination hasn't occured yet, add it to unique reads with UMI moved to read name, update occured reads
			if (seq1, seq2) not in occured:
				unique.append('%s;%s\n%s\n%s\n%s\n' % (info1, seq1[:args.umi_length], seq1[args.umi_length:], extra1, quality1[args.umi_length:]))
				occured.add((seq1, seq2))

	with gzip.open(args.output ,'wb') as out:
		out.write(''.join(unique))

def parseArguments():
	parser = argparse.ArgumentParser(prog="compress_UMI", description='Filters paired fastq file such that only a single instance of each read 1/read 2 sequence combination is kept. Moves UMI to read name.', usage='%(prog)s -n <UMI_length> -1 <read_1>.fastq.gz -2 <read_2>.fastq.gz -o <output>')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-n', '--UMI_length', type=int, required=True, help=' Length of the Unique Molecular Index', metavar='', dest='umi_length')
	required.add_argument('-1', '--read_1', required=True, help=' File containing read 1 (fastq.gz)', metavar='', dest='read_1')
	required.add_argument('-2', '--read_2', required=True, help=' File contianing read 2 (fastq.gz)', metavar='', dest='read_2')
	required.add_argument('-o', '--output', required=True, help=' Output file name', metavar='', dest='output')

	return parser.parse_args()

args = parseArguments()
main(args)
