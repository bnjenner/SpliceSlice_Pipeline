# SpliceSlice Analysis Pipeline

Pipeline for identification and statistical analysis of Branch Point and Polypyrimidine Tract sequences within introns. Method for extending BP_PPT.py to non-human organisms was modeled after this [publication](https://www.nature.com/articles/s42003-021-02725-7#Sec2), so extending this method to distantly related organisms should be approached with caution. Statistical analysis consists of a permutation test where the test statistic is the Kullback-Leibler (KL) divergence. Position weight matricies are also used to generate sequence logos for the identified sequences within each treatment.

Utilizes:
* [Bedtools](http://bedtools.readthedocs.org/) 
* [BP_PPT.py](https://github.com/zhqingit/BPP)
* [Scipy](https://docs.scipy.org/doc/scipy/index.html)
* [Logomaker](https://logomaker.readthedocs.io/en/latest/)
* Bash and Python Scripts


### Installation
1. Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for your System 

2. Clone this Repository
```
git clone https://github.com/bnjenner/SpliceSlice_Pipeline.git
```

3. Create Conda Environment
```
cd SpliceSlice_Pipeline
conda env create -f SpliceSlice.yaml
```

4. Activate Conda Environment
```
conda activate SpliceSlice
```

5. Add SpliceSlice.sh to your Path
Best to search how to do this on your specific system.

### Running SpliceSlice.sh
```
SpliceSlice.sh [-h] transcript_list1 transcript_list2 genome annotaion [ -o output_directory ]
```

Input files "transcript_list" should be a plain text file where each line contains one valid transcript ID and their quantification level (as an integer) separated by a tab.
