import sys

for line in open(sys.argv[1]):
	if line[0] != '#':
		cur = line.rstrip().split('\t')
		if cur[2] == 'CDS':
			print '%s\t%d\t%d' % (cur[0], int(cur[3])-1, int(cur[4])-1)