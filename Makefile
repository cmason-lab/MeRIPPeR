# MeRIPPeR-0.9.1
# Makefile
# Generated using MeRIPPeR->config.pl config.xml

# Sample Information
SAMPLES := Huh7_5 Huh7 TN_Huh7_5 TN_Huh7
Huh7_5 := Huh7_5_7_5
TN_Huh7_5 := TN_Huh7_5_7_5
Huh7_5_7_5_MERIP := Huh7_5_HCV_m6A
Huh7_5_7_5_CONTROL := Huh7_5_HCV_Flag
TN_Huh7_5_7_5_CONTROL := Huh7_5_HCV_m6A
TN_Huh7_5_7_5_MERIP := Huh7_5_HCV_Flag
Huh7 := Huh7_7
TN_Huh7 := TN_Huh7_7
Huh7_7_MERIP := Huh7_HCV_m6A
Huh7_7_CONTROL := Huh7_HCV_Flag
TN_Huh7_7_CONTROL := Huh7_HCV_m6A
TN_Huh7_7_MERIP := Huh7_HCV_Flag


# MeRIPPeR Configuration

REFSEQ_DIR := /home/yos2006/refseq
GENOME_FAI := hcv_genome/HCV.fasta.fai
BWA_INDEX := hcv_genome/bwa/HCV.fasta
GENOME_FA := hcv_genome/HCV.fasta
NUM_THREADS := 5
GENOME := HCV
ALPHA := 0.05
GSNAP_IIT := -s /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/GMAPDB/hg19/hg19.maps/hg19.refseq.splicesites.iit
WINDOW_SIZE := 25
GSNAP_GENOME_DIR := /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/GMAPDB
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
