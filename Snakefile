import os

# Build directory structure
for directory in ['full', 'downsample', 'downsample/fastq', 'logs', 'downsample/trim', 'logs/trim_report', 'logs/trim', 'downsample/align', 'logs/align', 'downsample/count', 'summary']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

strains = ['WT', 'prp2'] 
replicates = ['A', 'B', 'C']
sizes = [200000, 400000, 800000]
seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

rule all:
	input:
		'summary/downsampled_counts.txt',
		'summary/full_counts.txt'

###########################
# Full Dataset Processing #
###########################
rule trim_full:
	input:
		R1='raw/{strain}_{replicate}_R1.fastq.gz',
		R2='raw/{strain}_{replicate}_R2.fastq.gz'
	output:
		R1='full/trim/{strain}_{replicate}_R1.fastq.gz',
		R2='full/trim/{strain}_{replicate}_R2.fastq.gz',
		html='logs/trim_report/{strain}_{replicate}.html',
		json='logs/trim_report/{strain}_{replicate}.json'
	threads: 1
	log:
		'logs/trim/{strain}_{replicate}.log'
	params:
		'-w 1 --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT'
	shell:
		'fastp {params} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} 2> {log}'

rule align_full:
	input:
		R1='full/trim/{strain}_{replicate}_R1.fastq.gz',
		R2='full/trim/{strain}_{replicate}_R2.fastq.gz'
	output:
		bam = 'full/align/{strain}_{replicate}.bam'
	threads: 2
	log:
		'logs/align/{strain}_{replicate}.log'
	params:
		'--new-summary --max-intronlen 2000 --no-unal -x resources/hisat/sc_index_R64-2-1'
	shell:
		'hisat2 -p {threads} {params} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 5 | samtools sort -o {output.bam}'

rule feature_counts_full:
	input:
		'full/align/{strain}_{replicate}.bam'
	output:
		'full/count/{strain}_{replicate}.txt'
	threads: 1
	shell:
		'python2.7 scripts/feature_count.py --input {input} --output {output}'

rule combine_full:
	input:
		expand('full/count/{strain}_{replicate}.txt', strain=strains, replicate=replicates)
	output:
		'summary/full_counts.txt'
	threads: 1
	shell:
		'python scripts/combine_full.py --directory full/count/ --output {output}' 

###############################
# Downsampled File Processing #
###############################
rule downsample:
	input:
		R1 = 'raw/{strain}_{replicate}_R1.fastq.gz',
		R2 = 'raw/{strain}_{replicate}_R2.fastq.gz'
	output:
		expand('downsample/fastq/{{strain}}_{{replicate}}_{{seed}}_{size}_{read}.fastq.gz', size=sizes, read=['R1', 'R2'])
	threads: 1
	shell:
		'python scripts/downsample.py --input {wildcards.strain}_{wildcards.replicate} --in-dir raw/ --out-dir downsample/fastq --sizes %s --seed {wildcards.seed}' % (' '.join([str(size) for size in sizes])) 

rule trim_downsample:
	input:
		R1='downsample/fastq/{strain}_{replicate}_{seed}_{size}_R1.fastq.gz',
		R2='downsample/fastq/{strain}_{replicate}_{seed}_{size}_R2.fastq.gz'
	output:
		R1='downsample/trim/{strain}_{replicate}_{seed}_{size}_R1.fastq.gz',
		R2='downsample/trim/{strain}_{replicate}_{seed}_{size}_R2.fastq.gz',
		html='logs/trim_report/{strain}_{replicate}_{seed}_{size}.html',
		json='logs/trim_report/{strain}_{replicate}_{seed}_{size}.json'
	threads: 1
	log:
		'logs/trim/{strain}_{replicate}_{seed}_{size}.log'
	params:
		'-w 1 --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT'
	shell:
		'fastp {params} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} 2> {log}'

rule align_downsample:
	input:
		R1='downsample/trim/{strain}_{replicate}_{seed}_{size}_R1.fastq.gz',
		R2='downsample/trim/{strain}_{replicate}_{seed}_{size}_R2.fastq.gz'
	output:
		bam = 'downsample/align/{strain}_{replicate}_{seed}_{size}.bam'
	threads: 2
	log:
		'logs/align/{strain}_{replicate}_{seed}_{size}.log'
	params:
		'--new-summary --max-intronlen 2000 --no-unal -x resources/hisat/sc_index_R64-2-1'
	shell:
		'hisat2 -p {threads} {params} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 5 | samtools sort -o {output.bam}'

rule feature_count_downsample:
	input:
		'downsample/align/{strain}_{replicate}_{seed}_{size}.bam'
	output:
		'downsample/count/{strain}_{replicate}_{seed}_{size}.txt'
	threads: 1
	shell:
		'python2.7 scripts/feature_count.py --input {input} --output {output}'

rule combine_downsample:
	input:
		expand('downsample/count/{strain}_{replicate}_{seed}_{size}.txt', strain=strains, replicate=replicates, seed=seeds, size=sizes)
	output:
		'summary/downsampled_counts.txt'
	threads: 1
	shell:
		'python scripts/combine_downsample.py --directory downsample/count/ --output {output}' 
