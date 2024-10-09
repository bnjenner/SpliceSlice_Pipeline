import sys

# Listen, I know this isn't gonna be fast. 

with open(sys.argv[1], "r") as gi:

	gtf_lines = gi.readlines()
	gtf_dict = {}

	for i in range(len(gtf_lines)):

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



with open(sys.argv[2], "r") as fi:
	lines = fi.readlines()

for l in lines:

	if not l.startswith("PB_ID"):

		cols = l.split("\t")

		transcript_id = cols[3]
		pb_id = cols[0]

		if transcript_id == "novel":
			continue


		if gtf_dict[transcript_id]["strand"] == "+":

			# necessary copy, I am sorry
			introns = gtf_dict[transcript_id]["introns"][:-1]

			for x in range(0, len(introns)):

				# can we please stop with the different coordinate system
				#	or at least stop with open intervals, jc

				if (x % 2 == 0) or ((introns[x] - introns[x - 1]) < 300):
					continue

				print(gtf_dict[transcript_id]["chrom"],
					  introns[x - 1],
					  introns[x] - 1,
					 str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".full",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


				print(gtf_dict[transcript_id]["chrom"],
					  introns[x] - 1 - 252,
					  introns[x] - 1 - 1,
					  str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".end",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

				print(gtf_dict[transcript_id]["chrom"],
					  introns[x - 1],
					  introns[x - 1] + 50 ,
					 str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".5prime",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


				five_offset = -1 - 2 - 30 - 1 - 170 - 1
				three_offset = -1 - 2 - 30 - 1 # 0-based offset + dinucl offset + ppt window + end exclusion 

				print(gtf_dict[transcript_id]["chrom"],
					  introns[x] + five_offset,
					  introns[x] + three_offset,
					  str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".back",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


				five_offset = -1 - 2 - 30 - 1 # 0-based offset + dinucl offset + ppt window + end exclusion 
				three_offset = -1 - 2 

				print(gtf_dict[transcript_id]["chrom"], 
					  introns[x] + five_offset,
					  introns[x] + three_offset,
					  str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".ppt",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


		else:

			# if transcript_id == "ENSMUST00000162897.2":
			# 	print(gtf_dict[transcript_id])

			# necessary copy, I am sorryl
			introns = gtf_dict[transcript_id]["introns"][::-1]

			for x in range(0, len(introns)):

				if (x % 2 == 0) or ((introns[x + 1] - introns[x]) < 300):
					continue

				print(gtf_dict[transcript_id]["chrom"],
					  introns[x],
					  introns[x + 1] - 1,
					  str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".full",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

				print(gtf_dict[transcript_id]["chrom"],
					  introns[x] + 2,
					  introns[x] + 252 + 1,
					  str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".end",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

				print(gtf_dict[transcript_id]["chrom"],
					  introns[x + 1] - 1 - 50,
					  introns[x + 1] - 1,
					 str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".5prime",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

				

				five_offset = -1 + 3 + 30 + 1 # 0-based offset + dinucl offset + ppt window + end exlcusion
				three_offset = -1 + 3

				print(gtf_dict[transcript_id]["chrom"], 
					  introns[x] + three_offset,
					  introns[x] + five_offset,
					  str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".ppt",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


				five_offset = -1 + 3 + 30 + 2 + 170 # 0-based offset + dinucl offset + ppt window + end exlcusion + bps window
				three_offset = -1 + 3 + 30 + 1 # 0-based offset + dinucl offset + ppt window + end exlcusion

				print(gtf_dict[transcript_id]["chrom"], 
					  introns[x] + three_offset,
					  introns[x] + five_offset,
					  str(transcript_id) + "_" + pb_id + "_intron_" + str(x // 2) + ".back",
					  "0", gtf_dict[transcript_id]["strand"], sep = "\t")




# for transcript_id in gtf_dict.keys():

# 	if gtf_dict[transcript_id]["strand"] == "+":

# 		# necessary copy, I am sorry
# 		introns = gtf_dict[transcript_id]["introns"][:-1]

# 		for x in range(0, len(introns)):

# 			# can we please stop with the different coordinate system
# 			#	or at least stop with open intervals, jc

# 			if (x % 2 == 0) or ((introns[x] - introns[x - 1]) < 300):
# 				continue

# 			print(gtf_dict[transcript_id]["chrom"],
# 				  introns[x - 1],
# 				  introns[x] - 1,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".full",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

# 			print(gtf_dict[transcript_id]["chrom"],
# 				  introns[x] - 1 - 252,
# 				  introns[x] - 1 - 1,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".end",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


# 			five_offset = -1 - 2 - 30 - 1 - 170 - 1
# 			three_offset = -1 - 2 - 30 - 1 # 0-based offset + dinucl offset + ppt window + end exclusion 

# 			print(gtf_dict[transcript_id]["chrom"],
# 				  introns[x] + five_offset,
# 				  introns[x] + three_offset,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".back",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


# 			five_offset = -1 - 2 - 30 - 1 # 0-based offset + dinucl offset + ppt window + end exclusion 
# 			three_offset = -1 - 2 

# 			print(gtf_dict[transcript_id]["chrom"], 
# 				  introns[x] + five_offset,
# 				  introns[x] + three_offset,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".ppt",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


# 	else:

# 		# if transcript_id == "ENSMUST00000162897.2":
# 		# 	print(gtf_dict[transcript_id])

# 		# necessary copy, I am sorryl
# 		introns = gtf_dict[transcript_id]["introns"][::-1]

# 		for x in range(0, len(introns)):

# 			if (x % 2 == 0) or ((introns[x + 1] - introns[x]) < 300):
# 				continue

# 			print(gtf_dict[transcript_id]["chrom"],
# 				  introns[x],
# 				  introns[x + 1] - 1,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".full",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

# 			print(gtf_dict[transcript_id]["chrom"],
# 				  introns[x] + 2,
# 				  introns[x] + 252 + 1,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".end",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

			

# 			five_offset = -1 + 3 + 30 + 1 # 0-based offset + dinucl offset + ppt window + end exlcusion
# 			three_offset = -1 + 3

# 			print(gtf_dict[transcript_id]["chrom"], 
# 				  introns[x] + three_offset,
# 				  introns[x] + five_offset,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".ppt",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")


# 			five_offset = -1 + 3 + 30 + 2 + 170 # 0-based offset + dinucl offset + ppt window + end exlcusion + bps window
# 			three_offset = -1 + 3 + 30 + 1 # 0-based offset + dinucl offset + ppt window + end exlcusion

# 			print(gtf_dict[transcript_id]["chrom"], 
# 				  introns[x] + three_offset,
# 				  introns[x] + five_offset,
# 				  str(transcript_id) + "_intron_" + str(x // 2) + ".back",
# 				  "0", gtf_dict[transcript_id]["strand"], sep = "\t")

