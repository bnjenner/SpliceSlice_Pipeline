#!/bin/bash

# Function to display help message
usage() {
    echo "Usage: $0 [-h] transcript_list_file genome annotaion"
    echo
    echo "Positional arguments:"
    echo "  transcript_list_file    Path to file with list of transcripts of interest."
    echo "  genome                  Path to the genome (FASTA format)."
    echo "  annotation              Path to the genome annotation (GTF format)."
    echo "Optional arguments:"
    echo "  -o output_directory     Name of output directory."
    echo "  -h                      Show this help message and exit."
}

# Parse optional arguments
while getopts ":h" opt; do
    case ${opt} in
        h )
            usage
            exit 0
            ;;
        o )
			output=$OPTARG
			;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))


# Check for positional arguments
if [ $# -ne 3 ]; then
    echo "Error: Three positional arguments are required." 1>&2
    usage
    exit 1
fi

# Check for optional arguments
if [ -z "$output" ]; then
	output="SpliceSlice_Ouput"
fi

# Set Variables
transcript_file=$1
genome=$2
annotation=$3
transcript_file_prefix=`basename -s .txt $transcript_file`

# Slice Introns
echo "Slicing Introns..."
mkdir -p ${output}/00-IntronFiles
python3 scripts/intron_slicer.py -i $transcript_file \
						 -a $annotation --training \
						 > ${output}/00-IntronFiles/${transcript_file_prefix}.3prime.bed

# Get Sequence
bedtools getfasta -fi $genome \
          -bed ${output}/00-IntronFiles/${transcript_file_prefix}.3prime.bed \
          -s -name > ${output}/00-IntronFiles/${transcript_file_prefix}.3prime.fasta


# # Calculate Octanucleotide Frequences
# echo "Calculating Octanucleotide Frequences..."
# python3 scripts/get_freqy.py data/scPPT_human.txt \
#                 Splice_Slice/gencode.vM34.annotation.splice_slice.features.fasta \
#                 > Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt


# echo "Predicting Branch Point Sequences..."
# python scripts/BP_PPT.py -b Splice_Slice/pwmBP_human.txt \
#                  -p Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt \
#                  -i Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.intron_ends.fasta \
#                  > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.BP_predictions.txt


# # Set up Output Directory
# mkdir -p ${output}/00-IntronFiles
# mkdir -p ${output}/01-BP_Predictions
# mkdir -p ${output}/02-MotifAnalysis
