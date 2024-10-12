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

# Global
script_dir=$(dirname "$0")

##########################################################3
# Function 
slice() {

	local input=$1
	local ref=$2
	local gtf=$3
	local outpath=$4
	local prefix=$5
	local training=$6

	if [ $training == "true" ]; then
		python3 ${script_dir}/scripts/intron_slicer.py -i $input \
												 		-a $gtf \
												 		-o ${outpath}/${prefix}.introns.bed \
												 		--training \
														> ${outpath}/bp_training.bed

		bedtools getfasta -fi $ref \
					  	  -bed ${outpath}/bp_training.bed \
					  	  -s -name -fo ${outpath}/bp_training.fasta
	
	else 
		python3 ${script_dir}/scripts/intron_slicer.py -i $input \
								 		 			   -a $gtf \
								 		 			   -o ${outpath}/${prefix}.introns.bed
	
	fi

	# Get Sequences
	bedtools getfasta -fi $ref \
					  -bed <(grep ".target" ${outpath}/${prefix}.introns.bed) \
					  -s -name -fo ${outpath}/${prefix}.target.fasta

	# Get Sequences
	bedtools getfasta -fi $ref \
					  -bed <(grep ".ppt" ${outpath}/${prefix}.introns.bed) \
					  -s -name -fo ${outpath}/${prefix}.ppt.fasta	
}


freq() {

	local training=$1
	local outpath=$2

	python3 ${script_dir}/scripts/get_freqy.py ${script_dir}/data/scPPT_human.txt \
											   ${training} \
											   > ${outpath}/bp_training.octanucleotide_freqs.txt
}


predict() {

	local prefix=$1
	local inpath=$2
	local outpath=$3

	python ${script_dir}/scripts/BP_PPT.py -b ${script_dir}/data/pwmBP_human.txt \
                 			 			   -p ${inpath}/bp_training.octanucleotide_freqs.txt \
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


homer() {

	local group_1=$1
	local group_2=$2
	local outpath=$3

	echo "[     ${group_1} vs ${group_2}... ]" 

	# BP
	if [ ! -s ${outpath}/01-BP_Predictions/${group_1}.BP_predictions.fasta ]; then
		echo "[       No BP Predictions in ${outpath}/01-BP_Predictions/${group_1}.BP_predictions.fasta. ]"
	
	elif [ ! -s ${outpath}/01-BP_Predictions/${group_2}.BP_predictions.fasta ]; then
		echo "[       No BP Predictions in ${outpath}/01-BP_Predictions/${group_2}.BP_predictions.fasta. ]"
	
	else
		mkdir -p ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/BP/logs
		findMotifs.pl ${outpath}/01-BP_Predictions/${group_1}.BP_predictions.fasta \
					  fasta ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/BP -len 7 \
					  -fasta ${outpath}/01-BP_Predictions/${group_2}.BP_predictions.fasta \
					  1> ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/BP/logs/findMotifs.stdout \
					  2> ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/BP/logs/findMotifs.stderr
	fi

	# PPT
	if [ ! -s ${outpath}/00-IntronFiles/${group_1}.ppt.fasta ]; then
		echo "[       No PPT Predictions in ${outpath}/00-IntronFiles/${group_1}.ppt.fasta. ]"
	
	elif [ ! -s ${outpath}/00-IntronFiles/${group_2}.ppt.fasta ]; then
		echo "[       No PPT Predictions in ${outpath}/00-IntronFiles/${group_2}.ppt.fasta. ]"
	
	else 
		mkdir -p ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/PPT/logs
		findMotifs.pl ${outpath}/00-IntronFiles/${group_1}.ppt.fasta \
					  fasta ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/PPT -len 17 \
					  -fasta ${outpath}/00-IntronFiles/${group_2}.ppt.fasta \
					  1> ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/PPT/logs/findMotifs.stdout \
					  2> ${outpath}/02-Motif_Analysis/${group_1}_v_${group_2}/PPT/logs/findMotifs.stderr
	fi
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
slice $transcript_file_1 $genome $annotation ${output}/00-IntronFiles $prefix_1 "true"
slice $transcript_file_2 $genome $annotation ${output}/00-IntronFiles $prefix_2 "false"

	
# Calculate Octanucleotide Frequences
echo "[   Calculating Octanucleotide Frequences... ]"
freq ${output}/00-IntronFiles/bp_training.fasta ${output}/00-IntronFiles 


# Predict Branch Point Sequences
echo "[   Predicting Branch Point Sequences... ]"
mkdir -p ${output}/01-BP_Predictions
predict $prefix_1 ${output}/00-IntronFiles ${output}/01-BP_Predictions
predict $prefix_2 ${output}/00-IntronFiles ${output}/01-BP_Predictions


# Find Enriched Motifs
echo "[   Finding Enriched Sequence Motifs... ]"
homer $prefix_1 $prefix_2 $output
homer $prefix_2 $prefix_1 $output


end=`date +%s`
time=$((end-start))

echo "[ Finiahed in ${time} seoncds. ]"
echo "[ Pipeline Finished Successfully! ]"
