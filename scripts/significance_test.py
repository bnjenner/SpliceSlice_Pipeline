import argparse
import pandas as pd
import numpy as np
import itertools
import string
import logomaker
import math
import random
from collections import Counter
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



def generate_seq(pwm, length):
	seq = ''
	for i in range(length):
		probs = pwm.iloc[:, i]
		seq += random.choices(population=['A', 'C', 'G', 'T'], weights=probs, k=1)[0]

	return seq

def generate_3mers():
	combinations = itertools.product(["A", "C", "G", "T"], repeat=2)
	kmers = [''.join(combination) + "A" for combination in combinations]
	return kmers


def get_counts(file, kmers):

	with open(file) as fi:

		lines =	fi.readlines()

		sequences = []
		kmer_counts = pd.Series(1, index=kmers)

		n = len(lines[1].strip())
		data = np.ones((4, n)) # avoid weird zero division errors
		pwm = pd.DataFrame(data, index=['A', 'C', 'G', 'T'], columns=range(n))

		for seq in lines:
			if not seq.startswith('>'):
				sequences.append(seq.strip())
				if seq[3:6] in kmer_counts:
					kmer_counts.loc[seq[3:6]] += 1
				for i in range(n):
					pwm.loc[seq[i], i] += 1

		return pwm / sum(pwm), kmer_counts, sequences


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


def score_seq(pwm, seq):
	score = 0
	for i in range(len(seq)):
		score += math.log(pwm.loc[seq[i], i])

	return score


def KL_divergence(target, background):
	return stats.entropy(target, background)


def chisquare_test(observed, background):
	expected = (background / sum(background)) *  sum(observed)
	return stats.chisquare(observed, expected)


def mann_whitney_test(target_scores, generated_scores):
	return stats.mannwhitneyu(target_scores, background_scores)


##############################################################
if __name__ == "__main__":
	
	args = parse_arguments()

	kmers = generate_3mers()


	target_pwm, target_kmers, target_sequences = get_counts(args.target, kmers)
	background_pwm, background_kmers, _ = get_counts(args.background, kmers)

	con_target = consensus(target_pwm)
	con_background = consensus(background_pwm)

	make_logo(target_pwm.T, args.output, "Target")
	make_logo(background_pwm.T, args.output, "Background")

	# Monte Carlo Simulation
	num_seqs = 1000
	length_seq = background_pwm.shape[1]
	target_gen_seqs = [generate_seq(target_pwm, length_seq) for i in range(num_seqs)]
	background_gen_seqs = [generate_seq(background_pwm, length_seq) for i in range(num_seqs)]


	target_scores = [score_seq(background_pwm, seq) for seq in target_gen_seqs]
	background_scores = [score_seq(background_pwm, seq) for seq in background_gen_seqs]

	stat, p_value = mann_whitney_test(target_scores, background_scores)

	# kl = KL_divergence(target_kmers, background_kmers)
	# stat, p_value = chisquare_test(target_kmers, background_kmers)
	# stat, p_value = fisher_test(np.array(target_kmers), np.array(background_kmers))

	# print(stat, p_value)

	# with open(args.output + ".results.txt", "w") as fo:
	# 	fo.write("Target_Consensus\tBackground_Consensus\tKL\tX2_Stat\tP_Value\n")
	# 	fo.write(con_target + "\t" + con_background + "\t" + str(kl) + "\t" + str(stat) + "\t" + str(p_value) + "\n")
	# 