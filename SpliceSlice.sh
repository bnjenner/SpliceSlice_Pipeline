#!/bin/bash

##########################################################3
# Function to display help message
usage() {
    echo "Usage: $0 transcript_list1 transcript_list2 genome annotaion [ -o output_directory ]"
    echo
    echo "Positional arguments:"
    echo "  transcript_list1        Path to file with first list of transcripts."
    echo "  transcript_list2        Path to file with second list of transcripts."
    echo "  genome                  Path to the genome (FASTA format)."
    echo "  annotation              Path to the genome annotation (GTF format)."
    echo "Optional arguments:"
    echo "  -o output_directory     Name of output directory."
    echo "  -h                      Show this help message and exit."
}

# Parse optional arguments
while getopts ":o:h" opt; do
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
if [ $# -ne 4 ]; then
    echo "Error: Three positional arguments are required." 1>&2
    usage
    exit 1
fi

# Check for optional arguments
if [ -z "$output" ]; then
	output="SpliceSlice_Ouput"
fi


##########################################################3
# Function 
slice() {

	local input=$1
	local ref=$2
	local gtf=$3
	local outpath=$4
	local prefix=$5

	python3 scripts/intron_slicer.py -i $input \
							 		 -a $gtf --training \
									 > ${outpath}/${prefix}.introns.bed

	# Get Sequences
	bedtools getfasta -fi $ref \
					  -bed <(grep ".target" ${outpath}/${prefix}.introns.bed) \
					  -s -name -fo ${outpath}/${prefix}.target.fasta

	# Temporary File needs to be created because process substitution does 
	#	not work and I have no idea why
	grep -v ".target" ${outpath}/${prefix}.introns.bed \
			> ${outpath}/${prefix}.training.bed
	bedtools getfasta -fi $ref \
					  -bed ${outpath}/${prefix}.training.bed \
					  -s -name -fo ${outpath}/${prefix}.training.fasta
}


freq() {

	local prefix=$1
	local outpath=$2

	python3 scripts/get_freqy.py \
				${outpath}/${prefix}.training.fasta \
                > ${outpath}/${prefix}.training.octanucleotide_freqs.txt
}


predict() {

	local prefix=$1
	local inpath=$2
	local outpath=$3

	python scripts/BP_PPT.py -b data/pwmBP_human.txt \
                 			 -p ${inpath}/${prefix}.training.octanucleotide_freqs.txt \
                 			 -i ${inpath}/${prefix}.target.fasta \
                 			 > ${outpath}/${prefix}.BP_predictions.txt


	if [ -f ${outpath}/${prefix}.BP_predictions.fasta ]; then
		> ${outpath}/${prefix}.BP_predictions.fasta
	fi

	while IFS= read -r line; do
		if [[ $line != \#* ]]; then
			echo $line | cut -d ' ' -f 1 >> ${outpath}/${prefix}.BP_predictions.fasta
			echo $line | cut -d ' ' -f 2 >> ${outpath}/${prefix}.BP_predictions.fasta
		fi
	done < ${outpath}/${prefix}.BP_predictions.txt

}

##########################################################3
# Pipeline
echo "[ SpliceSlice Analysis Pipeline ]"

start=`date +%s`

# Set Variables
transcript_file_1=$1
transcript_file_2=$2
genome=$3
annotation=$4

echo "[     Group_1: ${transcript_file_1} ]"
echo "[     Group_2: ${transcript_file_2} ]"
echo "[     Genome: ${genome} ]"
echo "[     Annotation: ${annotation} ]"

prefix_1=`basename -s .txt $transcript_file_1`
prefix_2=`basename -s .txt $transcript_file_2`

# Slice Introns
echo "[   Slicing Sets of Introns... ]"
mkdir -p ${output}/00-IntronFiles
slice $transcript_file_1 $genome $annotation \
		${output}/00-IntronFiles $prefix_1
slice $transcript_file_2 $genome $annotation \
		${output}/00-IntronFiles $prefix_2

	
# Calculate Octanucleotide Frequences
echo "[   Calculating Octanucleotide Frequences... ]"
freq $prefix_1 ${output}/00-IntronFiles 
freq $prefix_2 ${output}/00-IntronFiles



# Predict Branch Point Sequences
echo "[   Predicting Branch Point Sequences... ]"
mkdir -p ${output}/01-BP_Predictions
predict $prefix_1 ${output}/00-IntronFiles ${output}/01-BP_Predictions
predict $prefix_2 ${output}/00-IntronFiles ${output}/01-BP_Predictions


# Find Enriched Motifs
echo "[   Finding Enriched Sequence Motifs... ]"
mkdir -p ${output}/02-Motif_Analysis/${prefix_1}_v_${prefix_2}/logs
mkdir -p ${output}/02-Motif_Analysis/${prefix_2}_v_${prefix_1}/logs

# Group 1 vs Group 2
findMotifs.pl ${output}/01-BP_Predictions/${prefix_1}.BP_predictions.fasta \
			fasta ${output}/02-Motif_Analysis/${prefix_1}_v_${prefix_2} -len 7 \
			-fasta ${output}/01-BP_Predictions/${prefix_2}.BP_predictions.fasta \
			1> ${output}/02-Motif_Analysis/${prefix_1}_v_${prefix_2}/logs/findMotifs.stdout \
			2> ${output}/02-Motif_Analysis/${prefix_1}_v_${prefix_2}/logs/findMotifs.stderr


# Group 2 vs Group 1
findMotifs.pl ${output}/01-BP_Predictions/${prefix_2}.BP_predictions.fasta \
			fasta ${output}/02-Motif_Analysis/${prefix_2}_v_${prefix_1} -len 7 \
			-fasta ${output}/01-BP_Predictions/${prefix_1}.BP_predictions.fasta \
			1> ${output}/02-Motif_Analysis/${prefix_2}_v_${prefix_1}/logs/findMotifs.stdout \
			2> ${output}/02-Motif_Analysis/${prefix_2}_v_${prefix_1}/logs/findMotifs.stderr


end=`date +%s`
time=$((end-start))

echo "[ Finiahed in ${time} seoncds. ]"
echo "[ Pipeline Finished Successfully! ]"
