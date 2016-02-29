# MeRIPPeR-0.9.1
# Makefile
# Generated using MeRIPPeR->config.pl config.xml

# Sample Information
SAMPLES := Sample2 Sample1 TN_Sample2 TN_Sample1
Sample2 := Sample2_Replicate1
TN_Sample2 := TN_Sample2_Replicate1
Sample2_Replicate1_MERIP := m6A_IP2
Sample2_Replicate1_CONTROL := input_control2
TN_Sample2_Replicate1_CONTROL := m6A_IP2
TN_Sample2_Replicate1_MERIP := input_control2
Sample1 := Sample1_Replicate1
TN_Sample1 := TN_Sample1_Replicate1
Sample1_Replicate1_MERIP := m6A_IP1
Sample1_Replicate1_CONTROL := input_control1
TN_Sample1_Replicate1_CONTROL := m6A_IP1
TN_Sample1_Replicate1_MERIP := input_control1


# MeRIPPeR Configuration

REFSEQ_DIR := REFSEQ_DIR
GENOME_FAI := hcv_genome/HCV.fasta.fai
BWA_INDEX := hcv_genome/bwa/HCV.fasta
GENOME_FA := hcv_genome/HCV.fasta
NUM_THREADS := 5
GENOME := HCV
ALPHA := 0.05
GSNAP_IIT := -s GSNAP_DIR/hg19/hg19.maps/hg19.refseq.splicesites.iit
WINDOW_SIZE := 25
GSNAP_GENOME_DIR := GSNAP_DIR
TOPHAT_GENOME_DIR := hcv_genome/bowtie2/hcv
STAR_REMOVE :=  --genomeDir hcv_genome/star/ --genomeLoad Remove --outFileNamePrefix star_genome_remove
STAR_INDEX := hcv_genome/star/
ALIGNER := star
MAPPING_QUALITY_THRESHOLD := 20
WINDOW_MAX_SIZE := 200
WINDOW_MIN_SIZE := 100
STAR_LOAD :=  --genomeDir hcv_genome/star/ --genomeLoad LoadAndExit --outFileNamePrefix star_genome_load
GENOME_SIZES := hcv_genome/hcv.sizes
GENOME_DIR := hcv_genome/

########## DO NOT EDIT BELOW THIS LINE ##########
ifndef MERIPPER_HOME
MERIPPER_HOME := .
endif

include $(MERIPPER_HOME)/makefiles/Config.mk
