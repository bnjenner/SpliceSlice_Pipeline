import sys

# parse classifications
with open(sys.argv[1], "r") as fi:

	class_dict = {}
	class_lines = fi.readlines()

	for i in range(1, len(class_lines)):

		cols = class_lines[i].split("\t")

		class_dict[cols[0]] = {"gene": cols[6],
							   "struct_cat": cols[5],
							   "transcript": cols[7],
							   "exons": cols[4],
							   "adj.pval": "NA",
							   "significant": 0}
		

# parse DE file
with open(sys.argv[2], "r") as fi:
	
	DE_transcripts = {}
	DE_lines = fi.readlines()

	for i in range(1, len(DE_lines)):

		cols = DE_lines[i].split("\t")
		adj_pval = float(cols[6])
		transcript_id = cols[0]

		if adj_pval <= 0.05:
			class_dict[transcript_id]["significant"] = 1
			class_dict[transcript_id]["adj.pval"] = adj_pval


# New Header
# PB_ID, Gene, Structual_Category, Transrcipt, Transcript_DE, Pval, Exons
print("PB_ID\tGene_ID\tStructural_Category\tTranscript_ID\tDE\tAdj_Pval\tExons")
for k, v in class_dict.items():

	pb_id = k
	gene = class_dict[k]["gene"]
	struct_cat = class_dict[k]["struct_cat"]
	transcript_id = class_dict[k]["transcript"]
	de = class_dict[k]["significant"]
	adj_pval = class_dict[k]["adj.pval"]
	exons = class_dict[k]["exons"]

	if de != 0:

		print(pb_id, gene, struct_cat, transcript_id, de, adj_pval, exons, sep = "\t")


