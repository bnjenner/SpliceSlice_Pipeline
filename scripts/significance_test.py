import argparse
import pandas as pd
import numpy as np
import logomaker
import random
import matplotlib.pyplot as plt
from scipy import stats

##############################################################
# Parse Commandline Arguments
def parse_arguments():

    parser = argparse.ArgumentParser(
        description="This is a sample help message for a script."
    )

    parser.add_argument(
        "-t",
        "--target",
        type=str,
        required=True,
        help="Path to file of target sequences.",
    )
    parser.add_argument(
        "-b",
        "--background",
        type=str,
        required=True,
        help="Path to file of target sequences.",
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Prefix for output file."
    )

    args = parser.parse_args()
    return args


def read_fasta(file):

    with open(file) as fi:
        lines = fi.readlines()
        sequences = []

        for seq in lines:
            if not seq.startswith(">"):
                sequences.append(seq.strip())

        return sequences


# Create Postion Weigt Matrix
def get_pwm(sequences):

    base_index = {"A": 0, "C": 1, "G": 2, "T": 3}

    n = len(sequences[0])
    pwm = np.ones(4 * n)  # pseudocount

    for seq in sequences:
        for i in range(n):
            pwm[(i * 4) + base_index[seq[i]]] += 1

    for i in range(0, n * 4, 4):
        pwm[i : i + 4] = pwm[i : i + 4] / sum(pwm[i : i + 4])

    return pwm


# Get "Consensus" Sequence
def consensus(pwm):

    nucleotides = ["A", "C", "G", "T"]
    max_indices = np.argmax(pwm, axis=0)
    consensus = "".join([nucleotides[idx] for idx in max_indices])

    return consensus

# Get Sequence Logo
def make_logo(pwm, file, sequence):

    logo = logomaker.Logo(pwm)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_ylabel("Position")
    logo.ax.set_ylabel("Frequency")
    logo.ax.set_title(file.split("/")[-1] + " " + sequence)

    plt.savefig(file + "." + sequence + ".png")

# Calculate KL Divergence
def KL_divergence(target, background):
    kl_sum = 0
    for i in range(len(target) // 4):
        kl_sum += stats.entropy(
            target[(i * 4) : (i * 4) + 5], background[(i * 4) : (i * 4) + 5]
        )

    return kl_sum

# Perform Permutation Test
def permutation_test(target_sequences, background_sequences, test_kl, iterations):

    counts = 0
    mean_sum = 0
    len_target = len(target_sequences)
    combined_sequences = target_sequences + background_sequences

    # Permutation
    for i in range(iterations):

        random.shuffle(combined_sequences)
        rand_target = combined_sequences[:len_target]
        rand_background = combined_sequences[len_target:]
        kl = KL_divergence(get_pwm(rand_target), get_pwm(rand_background))

        mean_sum += kl

        # This has never hit there is no way this is working this well?
        if test_kl <= kl:
            counts += 1

    return counts / iterations, mean_sum / iterations


##############################################################
if __name__ == "__main__":

    args = parse_arguments()

    target_sequences = read_fasta(args.target)
    background_sequences = read_fasta(args.background)

    target_pwm = get_pwm(target_sequences)
    background_pwm = get_pwm(background_sequences)
    test_kl = KL_divergence(target_pwm, background_pwm)

    p_value, sim_kl_mean = permutation_test(
        target_sequences, background_sequences, test_kl, 1000
    )

    seq_length = len(target_sequences[0])
    target_df = pd.DataFrame(
        target_pwm.reshape(seq_length, 4),
        columns=["A", "C", "G", "T"],
        index=range(len(target_sequences[0])),
    )
    background_df = pd.DataFrame(
        background_pwm.reshape(seq_length, 4),
        columns=["A", "C", "G", "T"],
        index=range(len(target_sequences[0])),
    )

    # Create Consensus Sequence
    con_target = consensus(target_df.T)
    con_background = consensus(background_df.T)

    # Create Logos
    make_logo(target_df, args.output, "Target")
    make_logo(background_df, args.output, "Background")

    with open(args.output + ".results.txt", "w") as fo:
        fo.write("Target_Consensus\tBackground_Consensus\tKL\tSim_KL_Mean\tP_Value\n")
        fo.write(
            con_target
            + "\t"
            + con_background
            + "\t"
            + str(test_kl)
            + "\t"
            + str(sim_kl_mean)
            + "\t"
            + str(p_value)
            + "\n"
        )
