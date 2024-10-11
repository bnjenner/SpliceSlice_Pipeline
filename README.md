# SpliceSlice Analysis Pipeline

## In Development

Utilizes:
* Bedtools (http://bedtools.readthedocs.org/) 
* BP_PPT.py (https://github.com/zhqingit/BPP)
* HOMER (http://homer.ucsd.edu/homer/)
* Bash and Python Scripts


### Installation
1. Install Conda for your System 

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

Input files "transcript_list" should be a plain text file where each line contains one valid transcript ID.
