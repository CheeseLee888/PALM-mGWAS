
################################# parameter settings below #################################
##################################### required below #######################################
# PALM method; 1 for 'PALM' or 2 for 'PALM-mbQTL'
PALMmethod=1

# input and output folders (Relative to WORK)
inputFolder=example/input
outputFolder=example/output

# input files (Relative to inputFolder)
genoFile=geno
abdFile=abd.txt
covFile=cov.txt

# chromosome number to analyze (for step2 only); set to 'NULL' for all chromosomes
chrom=1

# other control
useOriginalFormat=TRUE

##################################### optional below #######################################
# Default paths for other files, can be modified as needed
SIF=palmmbqtl.sif
outputSnpFile=info_snp.txt
outputFeatureFile=info_feature.txt
outputSeqDepthFile=info_seqdepth.txt
palm1_step1_prefix=palm1_step1_allpheno
palm2_step1_prefix=palm2_step1_allpheno
if [[ "${chrom}" == "NULL" ]]; then
  palm1_step2_prefix=palm1_step2_allchr
  palm2_step2_prefix=palm2_step2_allchr
else
  palm1_step2_prefix=palm1_step2_chr${chrom}
  palm2_step2_prefix=palm2_step2_chr${chrom}
fi




###################### export variables for scripts (do not modify) ########################
SIF=${WORK}/${SIF}
inputFolder=${WORK}/${inputFolder}
outputFolder=${WORK}/${outputFolder}
genoFile=${inputFolder}/${genoFile}
abdFile=${inputFolder}/${abdFile}
covFile=${inputFolder}/${covFile}
outputSnpFile=${outputFolder}/${outputSnpFile}
outputFeatureFile=${outputFolder}/${outputFeatureFile}
outputSeqDepthFile=${outputFolder}/${outputSeqDepthFile}
palm1_step1_prefix=${outputFolder}/${palm1_step1_prefix}
palm2_step1_prefix=${outputFolder}/${palm2_step1_prefix}
palm1_step2_prefix=${outputFolder}/${palm1_step2_prefix}
palm2_step2_prefix=${outputFolder}/${palm2_step2_prefix}

export SIF
export PALMmethod
export inputFolder
export outputFolder
export genoFile
export abdFile
export covFile
export chrom
export useOriginalFormat
export outputSnpFile
export outputFeatureFile
export outputSeqDepthFile
export palm1_step1_prefix
export palm1_step2_prefix
export palm2_step1_prefix
export palm2_step2_prefix