# Computational methodology for the detection of CRBN neo-substrate

This is the GitHub repository associated with my master's thesis "Computational methodology for the detection of CRBN neo-substrate". The purpose of this repository is to provide the code needed to replicate the results of this master's thesis. Moreover, the validation sets of the experiments are also included.

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#prerequisites">Prerequisites</a></li>
    <li><a href="#preprocessing">Preprocessing</a></li>
    <li><a href="#scripts">Scripts</a></li>
    <li><a href="#results">Results</a></li>
  </ol>
</details>

## Prerequisites

### Installation

1. Clone this repository.
2. Create the conda environment and install required packages.


### Conda environment

The `maria_enviroment.yml` file contains all the dependencies used during the development of this project. You can create the corresponding conda environment by running the following command:

```bash
conda env create -f maria_environment.yml
```

### Datasets

If you would like to replicate the results, you will need the [validation set](https://github.com/MariaSantamera00/Computational-methodology-for-the-detection-of-CRBN-neo-substrate/tree/master/validation_set) that was used in the paper. Refer to the paper for further explanation of the origin of this set.


### Other requirements

- Python 3.6 or later
- PyTorch 1.9.0 or later
- Numpy 1.19.5 or later
- Pandas 1.2.4 or later


## Preprocessing

The validation set contains both raw data and preprocessed data. In case you want to start from zero with the raw data of this set or other data of interest, use the `pdb_treatment.sh` script for data processing and adaptation.

## Scripts 

Each of the scripts contained in this repository has a first line of information explaining its purpose and use. Refer to Appendix B of the paper for a more in depth explanation of each script. 

The [bin folder](https://github.com/MariaSantamera00/Computational-methodology-for-the-detection-of-CRBN-neo-substrate/tree/master/bin) contains different Python scripts used by the main scripts of the pipeline. 

### Queue scripts 
Scripts whose name ends with "_q.sh" are designed to be sent to queue. 

## Results
The results obtained after developing the process using the [validation dataset](https://github.com/MariaSantamera00/Computational-methodology-for-the-detection-of-CRBN-neo-substrate/tree/master/validation_set) and applying this procedure to the computational models obtained with AlphaFold2 for all the proteins of the human proteome can be found in the [results folder](https://github.com/MariaSantamera00/Computational-methodology-for-the-detection-of-CRBN-neo-substrate/tree/master/results). 

