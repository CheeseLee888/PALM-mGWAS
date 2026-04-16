#!/usr/bin/env bash

# Repository root
WORK="/absolute/path/to/PALM-mGWAS"

# Common output root
outputFolder="${WORK}/example/output/study1"

# Step0 input/output
abdFile="${WORK}/example/input/study1/abd.txt"
covFile="${WORK}/example/input/study1/cov.txt"
genoFile="${WORK}/example/input/study1/geno"
abdAlignedFile="${outputFolder}/abd_aligned.txt"
covAlignedFile="${outputFolder}/cov_aligned.txt"
outputSeqDepthFile="${outputFolder}/info_seqdepth.txt"

# Optional Step0 / Step1 settings
covarColList="AGE,SEX"
depthCol="NULL"
depth_filter="0"
prev_filter="0.1"
outputFeatureFile="${outputFolder}/info_feature.txt"

# Step1 output
NULLmodelFile="${outputFolder}/step1_allpheno.rda"

# Step2.1 output prefix
# Single all-chromosome run:
#   palm1_step2_prefix="${outputFolder}/step2_allchr"
# Chromosome-split run:
palm1_step2_prefix="${outputFolder}/step2_chr1"

# Optional Step2.1 settings
chrom="1"
featureColList="NULL"
minMAF="0.05"
minMAC="5"
outputSnpFile="${outputFolder}/info_snp.txt"
useCluster="FALSE"

# Step2.2 merge
# If omitted, step2_2.sh falls back to palm1_step2_prefix as input and
# step2InputPrefix as output.
step2MergeInputPrefix="${outputFolder}/step2_chr1"
step2MergeOutputPrefix="${outputFolder}/step2_allchr"

# Step2.3 correction
step2InputPrefix="${outputFolder}/step2_allchr"
step2OutputPrefix="NULL"

# Step3 meta-analysis
studyDirFile="${WORK}/example/input/study_dirs.txt"
pattern="step2_allchr_.*[.]txt$"
metaDir="${WORK}/example/output/meta"
metaPrefix="step3_meta"

# Visualization
visPattern="step3_meta_.*\\.txt$"
plotDir="${WORK}/example/output/plot"
pheno="g_Bifidobacterium"
snp="chr1:14678191:A:C"
showMeta="TRUE"
showHet="TRUE"
plotMinP="NA"
pCut1="5e-8"
pCut2="5e-8"
plotWidth="NA"
plotHeight="NA"
