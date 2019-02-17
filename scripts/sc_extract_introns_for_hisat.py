import sys

for line in open(sys.argv[1]):
	if line[0] != '#':
		cur = line.rstrip().split('\t')
		if cur[2] == 'intron':
			print '%s\t%d\t%d\t%s' % (cur[0], int(cur[3])-1, int(cur[4])-1, cur[6])