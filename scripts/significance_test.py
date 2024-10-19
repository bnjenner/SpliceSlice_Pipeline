import argparse
import pandas as pd
import numpy as np
import logomaker
import matplotlib.pyplot as plt
from scipy.stats import chi2

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


def make_logo(pwm, file, sequence):
	
	logo = logomaker.Logo(pwm)
	logo.style_spines(visible=False)
	logo.style_spines(spines=['left', 'bottom'], visible=True)
	logo.ax.set_ylabel("Position")
	logo.ax.set_ylabel("Frequency")
	logo.ax.set_title(file.split("/")[-1] + " " + sequence)

	plt.savefig(file + "." + sequence + ".png")


def g_test(target, background):

	observed = np.array([target, background])

	total = observed.sum()
	row_totals = np.sum(observed, axis=1)
	col_totals = np.sum(observed, axis=0)

	expected = np.outer(row_totals, col_totals) / total
	g_stat = 2 * np.sum(observed * np.log(observed / expected))
	df = (observed.shape[0] - 1) * (observed.shape[1] - 1)
	p_value = chi2.cdf(g_stat, df)

	return g_stat, p_value



##############################################################
if __name__ == "__main__":
	
	args = parse_arguments()

	target_pwm = get_pwm(args.target)
	background_pwm = get_pwm(args.background)

	logo_target = make_logo(target_pwm.T, args.output, "Target")
	logo_target = make_logo(target_pwm.T, args.output, "Background")

	with open(args.output + ".g_test.txt", "w") as fo:

		fo.write("Position\tG-Statistic\tP-Value\n")
		for i in range(target_pwm.shape[1]):

			g_stat, pval = g_test(target_pwm[i], background_pwm[i])
			fo.write(str(i) + "\t" + str(g_stat) + "\t" + str(pval) + "\n")

