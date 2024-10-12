import sys

c_sum = 0
w_sum = 0

def normalized(weight, w_sum):
	if w_sum == 0:
		return 0
	return weight / w_sum

def weighted(cppt, cback, c_sum):
	if (cback * c_sum) == 0:
		return 0
	return (cppt * 2) / (cback * c_sum)


with open(sys.argv[1], "r") as fi:
	lines = fi.readlines()

	bp_dict = {}
	freq_dict = {}

	for l in lines:

		if not l.startswith("#"):

			cols = l.split("\t")
			bp_dict[cols[0]] = int(cols[2])
			freq_dict[cols[0]] = {"back": 0, "ppt": 0}



with open(sys.argv[2], "r") as fi:

	lines = fi.readlines()

	feature = ""

	for l in lines:

		if l.startswith('>'):
			feature = l.split(".")[3].split(":")[0]

		else:
			seq = l.strip()

			for i in range(0, len(seq) - 8) :

				octa = seq[i:i+8]
				freq_dict[octa][feature] += 1


c_sum = 0
for o, c in freq_dict.items():
	c_sum += c["ppt"]

w_sum = 0
for o, c in freq_dict.items():
	freq_dict[o]["weight"] = weighted(c["ppt"], c["back"], c_sum)
	w_sum += freq_dict[o]["weight"] 


print("#suseq\tn_ppt\tn_bp\tn_back\tscore")
for o, c in freq_dict.items():
	print(o, c["ppt"], bp_dict[o], c["back"], normalized(c["weight"], w_sum), sep = "\t")

