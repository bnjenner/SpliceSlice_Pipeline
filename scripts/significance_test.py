import argparse
import pandas as pd
import numpy as np
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


def get_pwm(file):

	with open(file) as fi:

		lines =	fi.readlines()

		n = len(lines[1].strip())
		data = np.ones((4, n)) # avoid weird zero division errors
		pwm = pd.DataFrame(data, index=['A', 'C', 'G', 'T'], columns=range(n))

		for seq in lines:
			if not seq.startswith('>'):
				for i in range(n):
					pwm.loc[seq[i], i] += 1

		return pwm 


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


def LRT(observed, expected):

	lrt = 2 * np.sum(observed * np.log(observed / expected))
	df = len(observed) - 1
	p_value = 1 - stats.chi2.cdf(lrt, df)

	return lrt, p_value



##############################################################
if __name__ == "__main__":
	
	args = parse_arguments()

	target_pwm = get_pwm(args.target)
	background_pwm = get_pwm(args.background)

	con_target = consensus(target_pwm)
	con_background = consensus(background_pwm)


	logo_target = make_logo(target_pwm.T, args.output, "Target")
	logo_target = make_logo(background_pwm.T, args.output, "Background")

	kl = KL_divergence(np.array(target_pwm).flatten(), np.array(background_pwm).flatten())
	lrt, p_value = LRT(np.array(target_pwm).flatten(), np.array(background_pwm).flatten())

	print("Target_Consensus\tBackground_Consensus\tKL\tLRT\tP-Value")
	print(con_target + "\t" + con_background + "\t" + str(kl) + "\t" + str(lrt) + "\t" + str(p_value))
	# with open(args.output + ".g_test.txt", "w") as fo:

	# 	fo.write("Position\tG-Statistic\tP-Value\n")
	# 	for i in range(target_pwm.shape[1]):

	# 		g_stat, pval = g_test(target_pwm[i], background_pwm[i])
	# 		fo.write(str(i) + "\t" + str(g_stat) + "\t" + str(pval) + "\n")

