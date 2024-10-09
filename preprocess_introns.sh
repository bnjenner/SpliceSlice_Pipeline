#!/bin/bash

echo "Getting Transcript Summary..."
python3 transcript_summary.py sqanti3/clustered.aligned.collapsed_classification.filtered_lite_classification.txt \
	2024_Janurary_Oqani_Mouse_Isoseq_R_analysis/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.txt \
	> 2024_Janurary_Oqani_Mouse_Isoseq_R_analysis/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.transcript_summary.txt

python3 transcript_summary.py sqanti3/clustered.aligned.collapsed_classification.filtered_lite_classification.txt \
        2024_Janurary_Oqani_Mouse_Isoseq_R_analysis/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.txt \
        > 2024_Janurary_Oqani_Mouse_Isoseq_R_analysis/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.transcript_summary.txt


echo "Getting Intron Bed..."
python3 splice_slicer.py References/gencode.vM34.annotation.gtf \
			2024_Janurary_Oqani_Mouse_Isoseq_R_analysis/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.transcript_summary.txt \
			> Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.introns.bed

python3 splice_slicer.py References/gencode.vM34.annotation.gtf \
                        2024_Janurary_Oqani_Mouse_Isoseq_R_analysis/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.transcript_summary.txt \
                        > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.introns.bed


grep ".end" Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.introns.bed \
	 > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.intron_ends.bed
grep ".ppt" Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.introns.bed \
         > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.ppt.bed
grep ".5prime" Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.introns.bed \
         > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.5prime.bed


grep ".end" Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.introns.bed \
         > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.intron_ends.bed
grep ".ppt" Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.introns.bed \
         > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.ppt.bed
grep ".5prime" Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.introns.bed \
         > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.5prime.bed



echo "Getting Intron Sequences"
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
          -bed Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.intron_ends.bed \
          -s -name > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.intron_ends.fasta
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
          -bed Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.ppt.bed \
          -s -name > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.ppt.fasta
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
          -bed Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.5prime.bed \
          -s -name > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.5prime.fasta


bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
          -bed Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.intron_ends.bed \
          -s -name > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.intron_ends.fasta
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
          -bed Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.ppt.bed \
          -s -name > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.ppt.fasta
bedtools getfasta -fi ../../References/GRCm39.primary_assembly.genome.fa \
          -bed Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.5prime.bed \
          -s -name > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.5prime.fasta



echo "Predicting Branch Point Sequences"
python BP_PPT.py -b Splice_Slice/pwmBP_human.txt \
                 -p Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt \
                 -i Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.intron_ends.fasta \
                 > Splice_Slice/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.BP_predictions.txt

python BP_PPT.py -b Splice_Slice/pwmBP_human.txt \
                 -p Splice_Slice/gencode.vM34.annotation.splice_slice.octanucleotide_freqs.txt \
                 -i Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.intron_ends.fasta \
                 > Splice_Slice/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.BP_predictions.txt

