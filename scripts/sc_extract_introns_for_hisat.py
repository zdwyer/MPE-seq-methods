import sys

def main(args):

	for line in open(sys.argv[1]):
		if line[0] != '#':
			cur = line.rstrip().split('\t')
			if cur[2] == 'intron':
				print '%s\t%d\t%d\t%s' % (cur[0], int(cur[3])-1, int(cur[4])-1, cur[6])

def parseArguments():
	parser = argparse.ArgumentParser(prog="sc_extract_introns_for_hisat.py", description='Extract the intron ranges from a feature file in a format suitable for hisat', usage='%(prog)s -f <feature file>.gff')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-f', '--feature_file', required=True, help=' Feature file', metavar='', dest='feature_file')
	return parser.parse_args()

args = parseArguments()
main(args)