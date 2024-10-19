import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

##############################################################
# Parse Commandline Arguments
def parse_arguments():

	parser = argparse.ArgumentParser(description='This is a sample help message for a script.')

	parser.add_argument('-t', '--target', type=str, required=True, help='Path to file of target sequences.')
	parser.add_argument('-b', '--background', type=str, required=True, help='Path to file of target sequences.')

	args = parser.parse_args()
	return args


def get_pwm(file):

	with open(file) as fi:

		lines =	fi.readlines()

		n = len(lines[1].strip())
		data = np.zeros((4, n))
		pwm = pd.DataFrame(data, index=['A', 'C', 'G', 'T'], columns=range(n))

		for seq in lines:
			if not seq.startswith('>'):
				for i in range(n):
					pwm.loc[seq[i], i] += 1

		col_sum = pwm.sum(axis=0)
		pwm = pwm / col_sum

		return pwm


def chi_sqr(target_df, background_df):

	chi2, p, dof, ex = chi2_contingency(target_df, background_df)
	return p


##############################################################
if __name__ == "__main__":
	
	args = parse_arguments()

	target_df = get_pwm(args.target)
	background_df = get_pwm(args.background)

	result = chi_sqr(target_df, background_df)
	print("Chi-Square P-Value: " + str(result))
