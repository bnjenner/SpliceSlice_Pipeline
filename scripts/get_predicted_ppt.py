import sys

with open(sys.argv[2], 'r') as fi:

	lines = fi.readlines()

	result = {}

	curr_header = ""

	for line in lines:

		if line.startswith('>'):
			curr_header = line.strip()
			result[curr_header] = ""

		else:

			if curr_header:
				result[curr_header] += line.strip()


with open (sys.argv[1], 'r') as fi:

	lines = fi.readlines()

	for line in lines:

		if not line.startswith('#'):

			cols = line.split('\t')
			fasta_id = cols[0]
			target_location = fasta_id.split(":")[-1].split('(')
			chrom = fasta_id.split(":")[-2]
			target_range = target_location[0]
			target_strand = target_location[-1][0]

			split_range = target_range.split('-')
			target_start = int(split_range[0])
			target_stop = int(split_range[-1])

			bp_A_pos = int(cols[2])

			if target_strand == '+':
				ppt_start = target_stop - (bp_A_pos + 2)
				ppt_stop = min(ppt_start + 20, target_stop + 10) + 1
			else:
				ppt_start = max(target_start - 10, target_start + (bp_A_pos - 2) - 20)
				ppt_stop = ppt_start + 20 + 1
				
			print(chrom, ppt_start, ppt_stop, fasta_id[1:], 0, target_strand, sep = "\t") 	


