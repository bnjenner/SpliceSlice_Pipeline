import sys, argparse

##############################################################
# Parse Commandline Arguments
def parse_arguments():

	parser = argparse.ArgumentParser(description='This is a sample help message for a script.')

	parser.add_argument('-i', '--input', type=str, required=True, help='Path to file with list of transcripts of interest.')
	parser.add_argument('-a', '--annotation', type=str, required=True, help='Path to annotation file (GTF format).')
	parser.add_argument('-b', '--bp-window', type=int, default=200, required=False, help='Length of window for Branch Point sequences.')
	parser.add_argument('-u', '--utr-window', type=int, default=150, required=False, help='Length of window for 5\' UTR.')
	parser.add_argument('-o', '--output', type=str, required=False, help='Path to transcripts of interest output file.')
	parser.add_argument('--training', action='store_true', help='Extract BPs, PPTs, and Background sequences for training.')

	args = parser.parse_args()
	return args


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

"""
BP Target zone is 10 to [W + 10] (W = window size) bps upstream of downstream exon start
187–200 nt and 3–16nt.
BACK           PPT

https://academic.oup.com/bioinformatics/article/33/20/3166/3870482
"""

def get_bed(input_file, bp_window, utr_window, output, training):

	with open(input_file, "r") as fi:
		lines = fi.readlines()
		transcripts_of_interest = {}

		for l in lines:
			if not l.startswith("#"): # just in case
				ID = l.split("\t")[0].strip()
				score = l.split("\t")[1].strip()
				transcripts_of_interest[ID] = score

	output_directory = {}

	for transcript_id in gtf_dict.keys():

		if len(gtf_dict[transcript_id]["introns"]) == 1:
			# print('NOTICE: ' + transcript_id + ' has no introns', file=sys.stderr)
			pass

		if gtf_dict[transcript_id]["strand"] == "+":

			# Forward Strand
			introns = gtf_dict[transcript_id]["introns"] # necessary copy, I am sorry

			for x in range(0, len(introns)):

				# can we please stop with the different coordinate system
				#	or at least stop with open intervals

				if (x % 2 == 0) or ((introns[x] - introns[x - 1]) < (bp_window + 100)):
					continue

				ppt = ""

				if training and ((introns[x + 1] - introns[x]) < 300):

					# Background Training
					print(gtf_dict[transcript_id]["chrom"],
						  introns[x] - 1 - (200), 
						  introns[x] - 1 - (187) + 1, #  bed files are half open
						  str(transcript_id) + ".intron_" + str(x // 2) + ".back",
						  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

					# PPT Training
					ppt = [gtf_dict[transcript_id]["chrom"],
						  str(introns[x] - 1 - (16) - 2),
						  str(introns[x] - 1 - 2 + 1),
						  str(transcript_id) + ".intron_" + str(x // 2) + ".ppt",
						  "0", gtf_dict[transcript_id]["strand"]]
					print("\t".join(ppt))

				# Potential Targets
				if transcript_id in transcripts_of_interest:
					output_directory[transcript_id] = (gtf_dict[transcript_id]["chrom"] + "\t" +
														 str(introns[x] - 1 - (10 + bp_window)) + "\t" +
														 str(introns[x] - 1 - 10 + 1) + "\t" +
														 str(transcript_id) + ".intron_" + str(x // 2) + ".target" + 
														 "." + transcripts_of_interest[transcript_id] + "\t" +
														 "0" + "\t" + gtf_dict[transcript_id]["strand"])
					if ppt != "":
						ppt[3] = str(transcript_id) + ".intron_" + str(x // 2) + ".ppt" +  "." + transcripts_of_interest[transcript_id]
						output_directory[transcript_id] += "\n" + "\t".join(ppt)


		else:

			# Reverse Strand
			introns = gtf_dict[transcript_id]["introns"][::-1]

			for x in range(0, len(introns)):

				if (x % 2 == 0) or ((introns[x + 1] - introns[x]) < (bp_window + 100)):
					continue

				ppt = ""

				if training and ((introns[x + 1] - introns[x]) < 300):

					# Background Training
					print(gtf_dict[transcript_id]["chrom"],
						  introns[x] - 1 + (187),
						  introns[x] - 1 + (200) + 1,
						  str(transcript_id) + ".intron_" + str(x // 2) + ".back",
						  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

					# PPT Training
					ppt = [gtf_dict[transcript_id]["chrom"],
						  str(introns[x] - 1 + (3)),
						  str(introns[x] - 1 + (16) + 1),
						  str(transcript_id) + ".intron_" + str(x // 2) + ".ppt",
						  "0", gtf_dict[transcript_id]["strand"]]
					print("\t".join(ppt))

				# Potential Targets
				if transcript_id in transcripts_of_interest:
					output_directory[transcript_id] = (gtf_dict[transcript_id]["chrom"] + "\t" +
														 str(introns[x] - 1 + 10) + "\t" +
														 str(introns[x] - 1 + (10 + bp_window) + 1)  + "\t" + 
														 str(transcript_id) + ".intron_" + str(x // 2) + ".target" + 
														 "." + transcripts_of_interest[transcript_id] + "\t" +
														 "0" + "\t" + gtf_dict[transcript_id]["strand"])
					
					if ppt != "":
						ppt[3] = str(transcript_id) + ".intron_" + str(x // 2) + ".ppt" +  "." + transcripts_of_interest[transcript_id]
						output_directory[transcript_id] += "\n" + "\t".join(ppt)



	with open(output, "w") as fo:
		for t in output_directory.values():
			fo.write(t + "\n")


##############################################################
if __name__ == "__main__":
	
	args = parse_arguments()

	gtf_dict = parse_annotation(args.annotation)
	get_bed(args.input, args.bp_window, args.utr_window, args.output, args.training)
