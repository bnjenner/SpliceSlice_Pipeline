import argparse
import pandas as pd
import numpy as np
import itertools
import string
import logomaker
import math
import matplotlib.pyplot as plt
from scipy import stats

##############################################################
# Parse Commandline Arguments
def parse_arguments():

	parser = argparse.ArgumentParser(description='This is a sample help message for a script.')

	parser.add_argument('-t', '--target', type=str, required=True, help='Path to file of target sequences.')
	parser.add_argument('-b', '--background', type=str, required=True, help='Path to file of target sequences.')
	parser.add_argument('-o', '--output', type=str, required=False, help='Prefix for output file.')

	args = parser.parse_args()
	return args


def generate_3mers():
	combinations = itertools.product(["A", "C", "G", "T"], repeat=2)
	kmers = [''.join(combination) + "A" for combination in combinations]
	return kmers


def get_counts(file, kmers):

	with open(file) as fi:

		lines =	fi.readlines()

		kmer_counts = pd.Series(1, index=kmers)

		n = len(lines[1].strip())
		data = np.ones((4, n)) # avoid weird zero division errors
		pwm = pd.DataFrame(data, index=['A', 'C', 'G', 'T'], columns=range(n))

		for seq in lines:
			if not seq.startswith('>'):

				if seq[3:6] in kmer_counts:
					kmer_counts.loc[seq[3:6]] += 1
				for i in range(n):
					pwm.loc[seq[i], i] += 1

		return pwm, kmer_counts


def consensus(pwm):

	nucleotides = ['A', 'C', 'G', 'T']
	max_indices = np.argmax(pwm, axis=0)
	consensus = ''.join([nucleotides[idx] for idx in max_indices])

	return consensus


def make_logo(pwm, file, sequence):
	
	logo = logomaker.Logo(pwm)
	logo.style_spines(visible=False)
	logo.style_spines(spines=['left', 'bottom'], visible=True)
	logo.ax.set_ylabel("Position")
	logo.ax.set_ylabel("Frequency")
	logo.ax.set_title(file.split("/")[-1] + " " + sequence)

	plt.savefig(file + "." + sequence + ".png")


def KL_divergence(target, background):
	return stats.entropy(target, background)


def chisquare_test(observed, background):
	expected = (background / sum(background)) *  sum(observed)
	return stats.chisquare(observed, expected)


##############################################################
if __name__ == "__main__":
	
	args = parse_arguments()

	kmers = generate_3mers()


	target_pwm, target_kmers = get_counts(args.target, kmers)
	background_pwm, background_kmers = get_counts(args.background, kmers)

	con_target = consensus(target_pwm)
	con_background = consensus(background_pwm)

	make_logo(target_pwm.T, args.output, "Target")
	make_logo(background_pwm.T, args.output, "Background")


	kl = KL_divergence(target_kmers, background_kmers)
	stat, p_value = chisquare_test(target_kmers, background_kmers)
	stat, p_value = fisher_test(np.array(target_kmers), np.array(background_kmers))

	print(stat, p_value)

	with open(args.output + ".results.txt", "w") as fo:
		fo.write("Target_Consensus\tBackground_Consensus\tKL\tX2_Stat\tP_Value\n")
		fo.write(con_target + "\t" + con_background + "\t" + str(kl) + "\t" + str(stat) + "\t" + str(p_value) + "\n")
	