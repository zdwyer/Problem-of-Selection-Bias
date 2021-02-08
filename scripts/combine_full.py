import os, argparse
from collections import defaultdict

def main(args):

	premature = defaultdict(int)
	mature = defaultdict(int)
	features = set() 
	samples = set()
	
	for infile in os.listdir(args.input_directory):
		sample = os.path.splitext(infile)[0]
		samples.add(sample)
		for line in open('%s%s' % (args.input_directory, infile)):
			cur = line.rstrip().split('\t')
			if cur[0] != "Intron":
				features.add(cur[0])
				mature[(sample, cur[0])] = int(cur[1])
				premature[(sample, cur[0])] = int(cur[2])

	with open(args.output_file, 'w') as out:
		out.write('Strain\tReplicate\tIntron\tMature\tPremature\n')
		for sample in sorted(samples):
			info = sample.split('_')
			for feature in sorted(features):
				out.write('%s\t%s\t%s\t%d\t%d\n' % (info[0], info[1], feature, mature[(sample, feature)], premature[(sample, feature)]))

def parseArguments():
	parser = argparse.ArgumentParser(prog="combine.py", description='Combines all files in a directory into single file.', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-d', '--directory', required=True, help=' Path to directory of files to combine.', metavar='', dest='input_directory')
	required.add_argument('-o', '--output', required=True, help=' Output file.', metavar='', dest='output_file')
	return parser.parse_args()

args = parseArguments()
main(args)