#!/bin/bash

# Get Branch Point and PPT ranges
echo "Slicing Genome..."
#python3 splice_slicer.py References/gencode.vM34.annotation.gtf > Splice_Slice/gencode.vM34.annotation.splice_slice.bed

grep ".full" Splice_Slice/gencode.vM34.annotation.splice_slice.bed > Splice_Slice/gencode.vM34.annotation.splice_slice.full_introns.bed
grep ".end" Splice_Slice/gencode.vM34.annotation.splice_slice.bed > Splice_Slice/gencode.vM34.annotation.splice_slice.intron_ends.bed
grep -v ".full\|.end" Splice_Slice/gencode.vM34.annotation.splice_slice.bed > Splice_Slice/gencode.vM34.annotation.splice_slice.features.bed


# Get Intron sequences
echo "Getting Intron Sequences..."
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
          -bed Splice_Slice/gencode.vM34.annotation.splice_slice.full_introns.bed \
          -s -name > Splice_Slice/gencode.vM34.annotation.splice_slice.full_introns.fasta

# Get Branch Point and PPT sequences
echo "Getting BP and PPT Sequences..."
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
	  -bed Splice_Slice/gencode.vM34.annotation.splice_slice.features.bed \
	  -s -name > Splice_Slice/gencode.vM34.annotation.splice_slice.features.fasta

# Get Intron End sequences
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
         -bed Splice_Slice/gencode.vM34.annotation.splice_slice.intron_ends.bed \
         -s -name > Splice_Slice/gencode.vM34.annotation.splice_slice.intron_ends.fasta


# Calculate weighted octanucleotide frequences
#echo "Calculating Weighted Frequences..."
python3 freqy.py ~/software/BPP/demo/scPPT_human.txt \
 		Splice_Slice/gencode.vM34.annotation.splice_slice.features.fasta \
		> Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt

# Predict Branch Point Sequences
echo "Performing BP Preidctions..."
conda activate python2.7
python BP_PPT.py -b Splice_Slice/pwmBP_human.txt \
		 -p Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt \
		 -i Splice_Slice/gencode.vM34.annotation.splice_slice.intron_ends.fasta \
		 > Splice_Slice/gencode.vM34.annotation.splice_slice.BP_predictions.txt
conda deactivate

python BP_PPT.py -b Splice_Slice/pwmBP_human.txt \
                 -p Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt \
                 -i Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.intron_ends.fasta \
                 > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.BP_predictions.txt

python BP_PPT.py -b Splice_Slice/pwmBP_human.txt \
                 -p Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt \
                 -i Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.intron_ends.fasta \
                 > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.BP_predictions.txt

echo "Done!"
