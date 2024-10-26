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


def read_fasta(file):

	with open(file) as fi:
		lines =	fi.readlines()
		sequences = []
		
		for seq in lines:
			if not seq.startswith('>'):
				sequences.append(seq.strip())
				
		return sequences


def get_pwm(sequences):

	n = len(sequences[0])
	data = np.ones((4, n)) # pseudocount
	pwm = pd.DataFrame(data, index=['A', 'C', 'G', 'T'], columns=range(n))

	for seq in sequences:
		for i in range(n):
			pwm.loc[seq[i], i] += 1

	return pwm / pwm.sum()


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
	kl_sum = 0
	for i in range(len(target.loc["A"])):
		kl_sum += stats.entropy(target.iloc[:,i], background.iloc[:,i])
	return kl_sum


def permutation_test(target_sequences, background_sequences, test_kl, iterations):

	counts = 0
	mean_sum = 0
	iterations = 100
	min_num = max(len(target_sequences), len(background_sequences))
	combined_sequences = target_sequences + background_sequences

	# Permutation
	for i in range(iterations):

		# Is this appropriate if list size is unequal?
		# Do I have enough data for this? How much sampling is approparite for this dataset?
		# sample witout replacement
		rand_target = random.sample(combined_sequences, min_num)
		rand_background = random.sample(combined_sequences, min_num)
		kl = KL_divergence(get_pwm(rand_target), get_pwm(rand_background))
		mean_sum += kl
		if test_kl <= kl:
			counts += 1

	return counts / iterations, mean_sum / iterations



def chisquare_test(observed, background):
	expected = (background / sum(background)) *  sum(observed)
	return stats.chisquare(observed, expected)


def mann_whitney_test(target_scores, generated_scores):
	return stats.mannwhitneyu(target_scores, background_scores)


##############################################################
if __name__ == "__main__":
	
	args = parse_arguments()


	target_sequences = read_fasta(args.target)
	background_sequences = read_fasta(args.background)

	target_pwm = get_pwm(target_sequences)
	background_pwm = get_pwm(background_sequences)
	test_kl = KL_divergence(target_pwm, background_pwm)

	p_value, sim_kl_mean = permutation_test(target_sequences, background_sequences, test_kl, 100)

	# Create Consensus Sequence
	con_target = consensus(target_pwm)
	con_background = consensus(background_pwm)

	# Create Logos
	make_logo(target_pwm.T, args.output, "Target")
	make_logo(background_pwm.T, args.output, "Background")

	print(p_value)

	with open(args.output + ".results.txt", "w") as fo:
		fo.write("Target_Consensus\tBackground_Consensus\tKL\tSim_KL_Mean\tP_Value\n")
		fo.write(con_target + "\t" + con_background + "\t" + str(test_kl) + "\t" + str(sim_kl_mean) + "\t" + str(p_value) + "\n")
	