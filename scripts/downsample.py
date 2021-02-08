import argparse, gzip, random
from itertools import izip_longest

def main(args):

	reads = []

	#Set seed of random number generator
	random.seed(args.seed)

	with gzip.open('%s/%s_R1.fastq.gz' % (args.indir, args.input)) as i1, gzip.open('%s/%s_R2.fastq.gz' % (args.indir, args.input)) as i2:
		for line1, line2 in izip_longest(i1, i2):
			
			#Read one entry from read 1
			info1 = line1.rstrip()
			seq1 = next(i1)
			extra1 = next(i1)
			quality1 = next(i1)

			umi = seq1[:7]

			#Read one entry from read 2
			info2 = line2.rstrip()
			seq2 = next(i2)
			extra2 = next(i2)
			quality2 = next(i2)

			# Give each read a random number from [0, 1)
			if len(seq1) > 40:
				reads.append((random.random() ,'%s;%s\n' % (info1, umi), seq1[7:], extra1, quality1[7:], '%s;%s\n' % (info2, umi), seq2, extra2, quality2))

	# Sort reads based on their random number
	reads.sort(key=lambda x: x[0])

	out1 = []
	out2 = []

	# Output first n reads (or all reads if n is greater than the number of reads)
	for i in range(0,min(len(reads), max(args.sizes))):
		out1.append(reads[i][1])
		out1.append(reads[i][2])
		out1.append(reads[i][3])
		out1.append(reads[i][4])

		out2.append(reads[i][5])
		out2.append(reads[i][6])
		out2.append(reads[i][7])
		out2.append(reads[i][8])

	for n in args.sizes:
		with gzip.open('%s/%s_%d_%d_R1.fastq.gz' % (args.outdir, args.input, args.seed, n), 'wb') as outfile:
			outfile.write(''.join(out1[:min(n, len(reads))*4]))
		with gzip.open('%s/%s_%d_%d_R2.fastq.gz' % (args.outdir, args.input, args.seed, n), 'wb') as outfile:
			outfile.write(''.join(out2[:min(n, len(reads))*4]))

def parseArguments():
	parser = argparse.ArgumentParser(prog="downsample.py", description='Downsamples paired-end fastq.gz file. Outputs .fastq.gz files with n reads.', usage='%(prog)s --input_1 <read_1> -2 <read_2> -n <numReadsToKeep> -o <output>')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', '--input', required=True, help=' Read 1 input file name', metavar='', dest='input')
	required.add_argument('-d', '--in-dir', required=True, help=' Input directory', metavar='', dest='indir')
	required.add_argument('-o', '--out-dir', required=True, help=' Out directory', metavar='', dest='outdir')
	required.add_argument('-n', '--sizes', nargs='+', type=int, required=True, help=' Number of reads to downsample to.', metavar='', dest='sizes')
	required.add_argument('-s', '--seed', type=int, required=True, help=' Seed for random number generator', metavar='', dest='seed')

	return parser.parse_args()

args = parseArguments()
main(args)
