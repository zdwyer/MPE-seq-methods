import argparse, gzip, random

def main(args):

	reads = []

	#Set seed of random number generator
	random.seed(args.seed)

	with gzip.open(args.input) as infile:
		for line1 in infile:
			
			#Read one entry from read 1
			info1 = line1
			seq1 = next(infile)
			extra1 = next(infile)
			quality1 = next(infile)

			# Give each read a random number from [0, 1)
			reads.append((random.random() ,info1, seq1, extra1, quality1))

	# Sort reads based on their random number
	reads.sort(key=lambda x: x[0])

	out = []

	# Output first n reads (or all reads if n is greater than the number of reads)
	for i in range(0,min(len(reads), args.numReads)):
		out.append(reads[i][1])
		out.append(reads[i][2])
		out.append(reads[i][3])
		out.append(reads[i][4])

	with gzip.open(args.output,'wb') as outfile:
		outfile.write(''.join(out))

def parseArguments():
	parser = argparse.ArgumentParser(prog="downsample_SE.py", description='Downsamples single-end fastq.gz file. Outputs .fastq.gz file with n reads.', usage='%(prog)s -i <input> -s <seed> -n <numReadsToKeep> -o <output>')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', '--input', required=True, help=' Input file name', metavar='', dest='input')
	required.add_argument('-o', '--output', required=True, help=' Output file name', metavar='', dest='output')
	required.add_argument('-n', '--numReads', type=int, required=True, help=' Number of reads to keep', metavar='', dest='numReads')
	required.add_argument('-s', '--seed', type=int, required=True, help=' Seed for random number generator', metavar='', dest='seed')

	return parser.parse_args()

args = parseArguments()
main(args)