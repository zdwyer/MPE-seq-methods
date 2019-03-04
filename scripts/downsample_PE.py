import argparse, gzip, random
from itertools import izip_longest

def main(args):

	reads = []

	#Set seed of random number generator
	random.seed(args.seed)

	with gzip.open(args.input_1) as i1, gzip.open(args.input_2) as i2:
		for line1, line2 in izip_longest(i1, i2):
			
			#Read one entry from read 1
			info1 = line1
			seq1 = next(i1)
			extra1 = next(i1)
			quality1 = next(i1)

			#Read one entry from read 2
			info2 = line2
			seq2 = next(i2)
			extra2 = next(i2)
			quality2 = next(i2)

			# Give each read a random number from [0, 1)
			reads.append((random.random() ,info1, seq1, extra1, quality1, info2, seq2, extra2, quality2))

	# Sort reads based on their random number
	reads.sort(key=lambda x: x[0])

	out1 = []
	out2 = []

	# Output first n reads (or all reads if n is greater than the number of reads)
	for i in range(0,min(len(reads), args.numReads)):
		out1.append(reads[i][1])
		out1.append(reads[i][2])
		out1.append(reads[i][3])
		out1.append(reads[i][4])

		out2.append(reads[i][5])
		out2.append(reads[i][6])
		out2.append(reads[i][7])
		out2.append(reads[i][8])


	with gzip.open(args.output_1,'wb') as outfile:
		outfile.write(''.join(out1))
	with gzip.open(args.output_2,'wb') as outfile:
		outfile.write(''.join(out2))

def parseArguments():
	parser = argparse.ArgumentParser(prog="downsample_PE.py", description='Downsamples paired-end fastq.gz file. Outputs .fastq.gz files with n reads.', usage='%(prog)s --input_1 <read_1> -2 <read_2> -n <numReadsToKeep> -o <output>')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', '--input_1', required=True, help=' Read 1 input file name', metavar='', dest='input_1')
	required.add_argument('-j', '--input_2', required=True, help=' Read 2 input file name', metavar='', dest='input_2')
	required.add_argument('-o', '--output_1', required=True, help=' Read 1 output file name', metavar='', dest='output_1')
	required.add_argument('-p', '--output_2', required=True, help=' Read_2 output file name', metavar='', dest='output_2')
	required.add_argument('-n', '--numReads', type=int, required=True, help=' Number of reads to keep', metavar='', dest='numReads')
	required.add_argument('-s', '--seed', type=int, required=True, help=' Seed for random number generator', metavar='', dest='seed')

	return parser.parse_args()

args = parseArguments()
main(args)