#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10-12
#SBATCH --mem=128000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --output=slurmout/homer_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmout/homer_%A.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnjenner@ucdavis.edu

hostname
klist

start=`date +%s`
echo $HOSTNAME

module load homer/4.10

export baseP=/share/biocore/projects/Rinaudo_P_UCSF/2024_Janurary_Oqani_Mouse_Isoseq
export cwd=$baseP/scripts
export seqP=$baseP/04-MotifAnalysis
export outP=$seqP/FB_exclude_IVF/BP_PPT


if [ ! -d $outP ]
then
  mkdir -p $outP
fi


findMotifs.pl ${seqP}/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.BP_predictions.fasta \
	fasta $outP -len 7 \
	-fasta ${seqP}/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.BP_predictions.fasta


findMotifs.pl ${seqP}/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.ppt.fasta \
        fasta $outP -len 20 \
        -fasta ${seqP}/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.ppt.fasta


export outP=$seqP/FB_exclude_IVF/5prime

if [ ! -d $outP ]
then
  mkdir -p $outP
fi


findMotifs.pl ${seqP}/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.5prime.fasta \
        fasta $outP \
        -fasta ${seqP}/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.5prime.fasta



export outP=$seqP/IVF_exclude_FB/BP_PPT

if [ ! -d $outP ]
then
  mkdir -p $outP
fi


findMotifs.pl ${seqP}/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.BP_predictions.fasta \
        fasta $outP -len 7 \
        -fasta ${seqP}/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.BP_predictions.fasta


findMotifs.pl ${seqP}/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.ppt.fasta \
        fasta $outP -len 20 \
        -fasta ${seqP}/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.ppt.fasta


export outP=$seqP/IVF_exclude_FB/5prime

if [ ! -d $outP ]
then
  mkdir -p $outP
fi



findMotifs.pl ${seqP}/FB_v_IVF_DE_analysis.IVF_exclude_FB.significant_DOWN.5prime.fasta \
        fasta $outP \
        -fasta ${seqP}/FB_v_IVF_DE_analysis.FB_exclude_IVF.significant_UP.5prime.fasta




end=`date +%s`
runtime=$((end-start))
echo $runtime


