import sys, argparse

##############################################################
# Parse Commandline Arguments
def parse_arguments():

	parser = argparse.ArgumentParser(description='This is a sample help message for a script.')

	parser.add_argument('-i', '--input', type=str, required=True, help='Path to file with list of transcripts of interest.')
	parser.add_argument('-a', '--annotation', type=str, required=True, help='Path to annotation file (GTF format).')
	parser.add_argument('-b', '--bp-window', type=int, default=100, required=False, help='Length of window for Branch Point sequences.')
	parser.add_argument('-u', '--utr-window', type=int, default=150, required=False, help='Length of window for 5\' UTR.')

	args = parser.parse_args()
	return args


'''
	Only consider introns that are at least 100 bps larger than the basepair window itself
	Locations:
		- BP: begin at positions 10-100 () bps upstream of downstream exon
			https://www.nature.com/articles/s42003-023-05513-7
		- PPT: 

'''

##############################################################
# Parse Annotation File
def parse_annotation(annotation_file):
	
	with open(annotation_file, "r") as gi:

		gtf_lines = gi.readlines()
		gtf_dict = {}

		for i in range(len(gtf_lines)):

			transcript_id = ""

			if not gtf_lines[i].startswith("#"):

				cols = gtf_lines[i].split("\t")
				anno = cols[-1].split("\"")
				chrom = cols[0]
				strand = cols[6]
				gene_id = anno[1]
				transcript_id = anno[3]

				if cols[2] == "transcript":
					gtf_dict[transcript_id] = {"chrom": chrom,
											   "strand": strand,
											   "introns": []}

				elif cols[2] == "exon":

					if next((a for a in anno if "exon_number 1;" in a), False):

						if strand == "+":
							gtf_dict[transcript_id]["introns"].append(int(cols[4]))
						else:
							gtf_dict[transcript_id]["introns"].append(int(cols[3]))
					
					else:

						if strand == "+":
							gtf_dict[transcript_id]["introns"].extend([int(cols[3]), int(cols[4])])
						else:
							gtf_dict[transcript_id]["introns"].extend([int(cols[4]), int(cols[3])])


	return gtf_dict


##############################################################
# Generate Intron Bed File
def get_bed(input_file, bp_window, utr_window):

	with open(input_file, "r") as fi:
		lines = fi.readlines()

	for l in lines:

		if not l.startswith("#"): # just in case
			transcript_id = l.strip()

			if transcript_id in gtf_dict:

				if len(gtf_dict[transcript_id]["introns"]) == 1:
					print('NOTICE: ' + transcript_id + ' has no introns', file=sys.stderr)

				if gtf_dict[transcript_id]["strand"] == "+":
					introns = gtf_dict[transcript_id]["introns"] # necessary copy, I am sorry

					for x in range(0, len(introns)):

						# can we please stop with the different coordinate system
						#	or at least stop with open intervals, jc

						if (x % 2 == 0) or ((introns[x] - introns[x - 1]) < (bp_window + 100)):
							continue

						print(gtf_dict[transcript_id]["chrom"],
							  introns[x] - 1 - (bp_window) - 2, 
							  introns[x] - 1 - 2 + 1, #  bed files are half open
							  str(transcript_id) + "_intron_" + str(x // 2) + ".end",
							  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

				else:
					introns = gtf_dict[transcript_id]["introns"][::-1]
					
					for x in range(0, len(introns)):

						if (x % 2 == 0) or ((introns[x + 1] - introns[x]) < (bp_window + 100)):
							continue

						print(gtf_dict[transcript_id]["chrom"],
							  introns[x] + 2,
							  introns[x] + (bp_window) + 2 + 1, # bed files are half open
							  str(transcript_id) + "_intron_" + str(x // 2) + ".end",
							  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


if __name__ == "__main__":
	
	args = parse_arguments()
	print(f"Input file: {args.input}")
	print(f"output file: {args.annotation}")

	gtf_dict = parse_annotation(args.annotation)
	
	get_bed(args.input, args.bp_window, args.utr_window)
