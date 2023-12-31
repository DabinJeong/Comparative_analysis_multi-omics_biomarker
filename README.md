# Table of contents
* [Project description](#Project-description)
* [Setup](#setup)
* [Run](#run)

# Project description
Build reproduciple pipleline of comparative analysis of <a href="https://github.com/DabinJeong/Multi-omics_biomarker"><u><b>multi-omics biomarker discovery project</b></u></a>.
Implemented the following tools in <a href="https://www.nextflow.io">nextflow</a>, which is a workflow manager.
- <a href="https://www.embopress.org/doi/pdf/10.15252/msb.20178124"> MOFA</a> (2018, Mol. Sys. Biol.)  <br>
  MOFA is originally designed for *unsupervised* integration of multi-omics data. As biomarker discovery problem is feature selection probelm in *supervised* task, we fed a logistic regression model with patient representation generated by MOFA to measure the predictive power of biomarkers.
- <a href="https://academic.oup.com/bioinformatics/article/37/16/2405/6129092"> iDRW</a> (2021, Bioinformatics) <br>
  iDRW is *supervised* multi-omics integration method that represents a patient with a pathway activity vector. To avoid information leakage, we modified iDRW (Code available at ./iDRW/R/get.iDRWP.R)and additionally implemented iDRW_test (Code available at ./iDRW/R/get.idRWP_test.R).
- <a href="https://pubmed.ncbi.nlm.nih.gov/30657866/"> DIABLO</a> (2019, Bioinformatics) <br>
  Multi-omics biomarker discovery tool based on canonical correlation analysis.
- <a href="https://www.nature.com/articles/s41467-021-23774-w"> MOGONET</a> (2021, Nat. Comm.) <br> 
  Multi-omics biomarker discovery tool based on graph neural network on the basis of multi-omics patient similarity network.


# Setup
## Build docker image
~~~
# Pull base image from docker hub
docker pull dabinjeong/cuda:10.1-cudnn7-devel-ubuntu18.04

# Build docker image
docker build --tag biomarker_comparison:0.1.1 .
~~~
## Install workflow manager: Nextflow
~~~
conda create -n biomarker python=3.9
conda activate biomarker
conda install -c bioconda nextflow=21.04.0
~~~

# Run
~~~
nextflow run comparison.nf -c comparison.config -with-docker biomarker_comparison:0.1.1
~~~
