#####################
# MeRIPPeR 			#
# Config.mk			#
# Yogesh Saletore 	#
#####################

######################################

# REFSEQ ANNOTATION
MAIN_DIR := $(ALIGNER)
ALIGNMENT_DIR := $(ALIGNER)/01-alignment
PEAKS_DIR := $(ALIGNER)/02-peaks
ANNOTATION := $(ALIGNER)/03-annotation

SHELL := /bin/bash
TMP_DIR := tmp
REFSEQ_EXONS   := $(REFSEQ_DIR)/$(GENOME).exons.bed
REFSEQ_GENES   := $(REFSEQ_DIR)/$(GENOME).genes.bed
REFSEQ_GTF     := $(REFSEQ_DIR)/$(GENOME).gtf
REFSEQ_REFGENE := $(REFSEQ_DIR)/$(GENOME).txt

PEAKS := $(addsuffix .$(ALIGNER).peaks.bed, $(SAMPLES))
PEAKS_SPLIT := $(addsuffix .$(ALIGNER).peaks-split.bed, $(SAMPLES))
STAR_REMOVE_TARGET :=

ifeq ($(ALIGNER), star)
	STAR_REMOVE_TARGET := star.remove
endif

define R2S
$($(strip $(1))_MERIP) $($(strip $(1))_CONTROL)
endef

define S2R
$(foreach S, $($(strip $(1))) , $(call R2S, $(S)))
endef
DIRNAME := $(shell echo $(CURDIR) | sed 's~.*/~~' )
FILES := $(foreach S, $(SAMPLES), $(call S2R, $(S)))
CUFFLINKS_TARGETS := $(addsuffix .$(ALIGNER).cufflinks.transcripts.gtf, $(FILES)) $(addsuffix .$(ALIGNER).cufflinks.isoforms.fpkm_tracking, $(FILES)) $(addsuffix .$(ALIGNER).cufflinks.genes.fpkm_tracking, $(FILES))
REFSEQ_EXONS_TARGETS := $(addsuffix .$(ALIGNER).refseq.exons.bed, $(FILES))
#REFSEQ_GENES_TARGETS := $(addsuffix .$(ALIGNER).refseq.genes.bed, $(FILES))
REFSEQ_GENES_RPKM_TARGETS := $(addsuffix .$(ALIGNER).refseq.genes.rpkm.bed, $(FILES))
BIGWIG_TARGETS := $(addsuffix .$(ALIGNER).bw, $(FILES))
WC_TARGETS := $(addsuffix .$(ALIGNER).filt.bam.rc.txt, $(FILES)) $(addsuffix .$(ALIGNER).filt.bam.frac.txt, $(FILES))  $(addsuffix .fastq.gz.rc.txt, $(FILES)) $(addsuffix .$(ALIGNER).bam.frac.txt, $(FILES))
BAM_TARGETS := $(addsuffix .$(ALIGNER).filt.bam, $(FILES))
TPP_TARGETS := $(addsuffix .tpp.n0.pdf, $(PEAKS_SPLIT)) $(addsuffix .$(ALIGNER).tpp.n0.pdf, $(SAMPLES)) $(DIRNAME).$(ALIGNER).tpp.n0.pdf $(addsuffix .tpp.n1.pdf, $(PEAKS_SPLIT)) $(addsuffix .$(ALIGNER).tpp.n1.pdf, $(SAMPLES)) $(DIRNAME).$(ALIGNER).tpp.n1.pdf $(addsuffix .tpp.n2.pdf, $(PEAKS_SPLIT)) $(addsuffix .$(ALIGNER).tpp.n2.pdf, $(SAMPLES)) $(DIRNAME).$(ALIGNER).tpp.n2.pdf

PEAKS_COUNTS_DIR := $(TMP_DIR)/counts
PEAKS_COUNTS := $(addsuffix .counts.txt, $(PEAKS_SPLIT)) $(addsuffix .counts.txt, $(PEAKS))

# $(addsuffix .$(ALIGNER).001-raw-windows.bed, $(1)): $(addsuffix .$(ALIGNER).filt.bam, $($(strip $(1))_MERIP) $($(strip $(1))_CONTROL) )
define REPLICATE_TO_SEQ
$(addsuffix .$(ALIGNER).005-windows-merged-filtered.bed, $(1)): $(addsuffix .$(ALIGNER).filt.bam, $($(strip $(1))_MERIP) $($(strip $(1))_CONTROL) )
endef

define SAMPLE_TO_REPLICATE
$(addsuffix .$(ALIGNER).peaks-compl.bed, $(1)): $(addsuffix .$(ALIGNER).006-windows-merged-filtered-compl.bed, $($(strip $(1))) )
$(foreach S, $($(strip $(1))) , $(eval $(call REPLICATE_TO_SEQ, $(S))))
endef

### TOOLS ###
CAT := cat
ZCAT := zcat

# ALIGNERS
BWA := bwa
BWA_ALN_OPTIONS := aln -t $(NUM_THREADS)
BWA_SAMSE_OPTIONS := samse -n 1

PADJUST := --p.adjust BenjaminiHochberg

GSNAP := gsnap
GSNAP_OPTIONS := --gunzip -D $(GSNAP_GENOME_DIR) -d $(GENOME) -t $(NUM_THREADS) -A sam --novelsplicing=1 -m 4 -B 5 $(GSNAP_IIT)

TOPHAT := tophat2 #--GTF  $(REFSEQ_GTF)
TOPHAT_OPTIONS ?= --output-dir $*-tophat_out --max-multihits 1 --GTF /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/ANNOTATION/hg19/refseq/2015-02-11/hg19.refflat.gtf --num-threads $(NUM_THREADS) $(TOPHAT_GENOME_DIR) $<

ifndef STAR_INDEX
STAR_INDEX := /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/STAR/STAR_2.3.0e/$(GENOME)
endif
STAR := STAR
STAR_LOAD := $(STAR_LOAD) #--genomeDir /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/STAR/STAR_2.3.0e/$(GENOME) --genomeLoad LoadAndExit --outFileNamePrefix star_genome_load
STAR_REMOVE := $(STAR_REMOVE) #--genomeDir /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/STAR/STAR_2.3.0e/$(GENOME) --genomeLoad Remove --outFileNamePrefix star_genome_remove
STAR_MAPPING_OPTIONS ?= --genomeDir $(STAR_INDEX) --runThreadN $(NUM_THREADS) --readFilesIn $< --outStd SAM --outFilterMultimapNmax 10 --readFilesCommand zcat --outFileNamePrefix $@. --genomeLoad LoadAndKeep --outSAMstrandField intronMotif 

CUFFLINKS := cufflinks
CUFFLINKS_OPTIONS ?= --no-faux-reads --num-threads $(NUM_THREADS) --GTF $(REFSEQ_DIR)/$(GENOME).gtf --output-dir $*-cufflinks $^

SAMTOOLS := samtools
BAMTOBED := bamToBed
COVERAGEBED := coverageBed

BEDTOOLS := bedtools
BEDTOOLS_MULTICOV := $(BEDTOOLS) multicov 
BEDTOOLS_MULTIINTER := $(BEDTOOLS) multiinter

ifdef JAVA_HOME
JAVA := $(JAVA_HOME)/bin/java
else
JAVA := java
endif
JAVA_OPTS := -Dfile.encoding=UTF-8
MERIPPER_CLASSPATH := -cp "$(MERIPPER_HOME)/jars/MeRIPPeR.jar:$(MERIPPER_HOME)/jars/commons-logging-1.1.1.jar:$(MERIPPER_HOME)/jars/picard-1.68.jar:$(MERIPPER_HOME)/jars/sam-1.68.jar:$(MERIPPER_HOME)/jars/commons-cli-1.2.jar"

MERIPPER_WINDOWCOUNTER ?= $(JAVA) $(JAVA_OPTS) -Xmx2g -Xms512m $(MERIPPER_CLASSPATH) edu.cornell.med.icb.masonlab.meripper.peakfinder.WindowCounter
MERIPPER_WINDOWCOUNTER_OPTIONS ?= --merip $(firstword $^) --control $(lastword $^) --output $@ --format sam --genome-sizes $(GENOME_SIZES) --num-threads $(NUM_THREADS) --window-size $(WINDOW_SIZE)

MERIPPER_WINDOW_PVALUE_FILTER ?= $(JAVA) $(JAVA_OPTS) -Xmx2g -Xms2g $(MERIPPER_CLASSPATH) edu.cornell.med.icb.masonlab.meripper.peakfinder.PValueWindowFilter
MERIPPER_WINDOW_PVALUE_FILTER_OPTIONS ?= --input $< --output $@ --alpha $(ALPHA) $(PADJUST) --genome-sizes $(GENOME_SIZES) --window-size $(WINDOW_SIZE)

MERIPPER_WINDOW_SIZE_FILTER ?= $(JAVA) $(JAVA_OPTS) -Xmx2g -Xms2g $(MERIPPER_CLASSPATH) edu.cornell.med.icb.masonlab.meripper.peakfinder.WindowSizeFilter
MERIPPER_WINDOW_SIZE_FILTER_OPTIONS ?= --input $< --output $@ --min $(WINDOW_MIN_SIZE)

MERIPPER_WINDOWSPLITTER ?= $(JAVA) $(JAVA_OPTS) -Xmx2g -Xms512m $(MERIPPER_CLASSPATH) edu.cornell.med.icb.masonlab.meripper.peakfinder.WindowSplitter
MERIPPER_WINDOWSPLITTER_OPTIONS ?= --input $< --output $@ --max $(WINDOW_MAX_SIZE)

MERIPPER_TRANSCRIPTOME_PROFILE_PLOT := $(JAVA) $(JAVA_OPTS) -Xmx2g -Xms512m $(MERIPPER_CLASSPATH) edu.cornell.med.icb.masonlab.jenotator.plotting.transcriptome.TranscriptomeProfilePlot
MERIPPER_TRANSCRIPTOME_PROFILE_PLOT := $(JAVA) -Xmx5g -jar /home/yos2006/MeRIPPeR-tpp.jar
MERIPPER_TRANSCRIPTOME_PROFILE_PLOT_OPTIONS ?= --plot5UTR --plotCDS --plot3UTR --input $(shell echo $(strip $^) | sed 's/ /,/g') --legend $(shell echo $(strip $^) | sed 's/ /,/g') --input-type $(shell printf 'bed %.0s' {1..$(words $^)} | sed -e "s/^ *//" | sed -e "s/ *$$//" | sed 's/ /,/g') --transcripts $(REFSEQ_REFGENE) --genome-sizes $(GENOME_SIZES) --transcripts-type refgene --prefix $(shell echo $(addsuffix -, $(addprefix tmp/$@-, $(strip $^))) | sed 's/ /,/g') --output $@ --R-col blue,red,deeppink,green4,black,aquamarine,blueviolet,navy,green,orange,deepskyblue,darkred,magenta,slateblue4 --R-lty solid,solid,solid,solid,solid,dashed,dashed,dashed,dashed,dashed,dashed,dashed,dashed,dashed,dashed --R-lwd 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 --num-threads $(NUM_THREADS)

GENOME_COVERAGE_BED ?= genomeCoverageBed
GENOME_COVERAGE_BED_OPTIONS ?= -bga -split -ibam $^ -g $(GENOME_SIZES) | sort -k1,1 -k2,2n > $@

BEDGRAPH_TO_BIGWIG := bedGraphToBigWig
BEDGRAPH_TO_BIGWIG_OPTIONS ?= $^ $(GENOME_SIZES) $@

COMPLEMENT_BED := complementBed
COMPLEMENT_BED_OPTIONS ?= -g $(GENOME_SIZES)

SUBTRACT_BED := subtractBed
SUBTRACT_BED_OPTIONS ?= 

MERGE_BED := mergeBed

SORT := sort
SORT_OPTIONS := -k 1,1 -k 2,2n

PERL := perl
RPKM_PL := $(MERIPPER_HOME)/scripts/perl/rpkm.pl
FRAC_PL := $(MERIPPER_HOME)/scripts/perl/frac.pl

### TARGETS ###
# KEEP THIS BLANK SO INTERMEDIATES ARE KEPT!!!
.SECONDARY: 

.PHONY: all clean bigwig refseq cufflinks peaks peaks-annotate $(STAR_REMOVE_TARGET) wc

wc: $(WC_TARGETS)

default: peaks

peaks: $(PEAKS_SPLIT)

cufflinks: $(CUFFLINKS_TARGETS)

refseq: $(REFSEQ_EXONS_TARGETS) $(WC_TARGETS) $(REFSEQ_GENES_RPKM_TARGETS) $(ALIGNER).rpkm.txt $(ALIGNER).read-counts.txt $(ALIGNER).library-size.txt $(ALIGNER).mapped.txt

peaks-annotate: $(TPP_TARGETS)

all: peaks peaks-annotate $(PEAKS_COUNTS) $(WC_TARGETS)

%.bam.bai: %.bam
	$(SAMTOOLS) index $<

# BIGWIGS
bigwig: $(BIGWIG_TARGETS)

%.bw: %.bedgraph
	$(BEDGRAPH_TO_BIGWIG) $(BEDGRAPH_TO_BIGWIG_OPTIONS)

%.bedgraph: %.filt.bam
	$(GENOME_COVERAGE_BED) $(GENOME_COVERAGE_BED_OPTIONS)

# BWA ALIGNMENT
%.bwa.sai: %.fastq.gz
	$(BWA) $(BWA_ALN_OPTIONS) $(BWA_INDEX) $< > $@

%.bwa.bam: %.bwa.sai %.fastq.gz
	$(BWA) $(BWA_SAMSE_OPTIONS) $(BWA_INDEX) $^ | $(SAMTOOLS) view -bS - | $(SAMTOOLS) sort - $(basename $@)

# GSNAP ALIGNMENT
%.gsnap.bam: %.fastq.gz
	$(GSNAP) $(GSNAP_OPTIONS) $< | $(SAMTOOLS) view -bS - | $(SAMTOOLS) sort - $(basename $@)

# TOPHAT ALIGNMENT
%.tophat.bam: %-tophat_out/accepted_hits.bam
	$(SAMTOOLS) view $^ | awk '$$3 != "*"' | $(SAMTOOLS) view -bS -t $(GENOME_FAI) - | $(SAMTOOLS) sort - $(basename $@)

%-tophat_out/accepted_hits.bam: %.fastq.gz
	$(TOPHAT) $(TOPHAT_OPTIONS)

%.star.bam: %.fastq.gz star.loaded
	$(STAR) $(STAR_MAPPING_OPTIONS) | $(SAMTOOLS) view -bS -t $(GENOME_FAI) - | $(SAMTOOLS) sort - $(basename $@)

# hack because STAR MAPQ is based on multimappers.
# since we don't allow them, this is pointless
%.star.filt.bam: %.star.bam
	ln -sf $< $@ 
#	$(SAMTOOLS) view -b -q $(MAPPING_QUALITY_THRESHOLD) $< > $@
#	$(CAT) <($(SAMTOOLS) view -H $<) <($(SAMTOOLS) view $< | awk 'BEGIN {OFS = "\t"; FS = "\t"} {$$2 = !and($$2,0x0040); print $$0 }') | $(SAMTOOLS) view -Sbt $(GENOME_FAI) - > $@

star.loaded:
	$(STAR) $(STAR_LOAD)
	touch $@

star.remove: $(BAM_TARGETS)
	if [ -f star.loaded ]; then $(STAR) $(STAR_REMOVE); rm star.loaded; fi 

# ALIGNMENT QUALITY FILTERING
%.filt.bam: %.bam
	$(SAMTOOLS) view -b -q $(MAPPING_QUALITY_THRESHOLD) $< | $(SAMTOOLS) sort - $(basename $@)

#
# PEAK FINDING
#
#
# PER REPLICATE PER SAMPLE

#%.001-raw-windows.bed:
#	$(MERIPPER_WINDOWCOUNTER) $(MERIPPER_WINDOWCOUNTER_OPTIONS)

#%.002-filtered-windows.bed: %.001-raw-windows.bed
#	$(MERIPPER_WINDOW_PVALUE_FILTER) $(MERIPPER_WINDOW_PVALUE_FILTER_OPTIONS)

#%.003-filtered-windows-sorted.bed: %.002-filtered-windows.bed
#	$(SORT) $(SORT_OPTIONS) $^ > $@

#%.004-windows-merged.bed: %.003-filtered-windows-sorted.bed
#	$(CAT) $^ | $(MERGE_BED) > $@

#%.005-windows-merged-filtered.bed: %.004-windows-merged.bed
%.005-windows-merged-filtered.bed:
	java -Xmx10g -jar jars/MeRIPPeR.jar --merip $(firstword $^) --control $(lastword $^) --output $@ --genome-sizes $(GENOME_SIZES) --num-threads $(NUM_THREADS) --window-size $(WINDOW_SIZE) --alpha $(ALPHA) $(PADJUST) --min-window $(WINDOW_MIN_SIZE) --step-size 25 --genes /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/ANNOTATION/hg19/refseq/2015-02-11/hg19.refflat.bed
#	$(MERIPPER_WINDOW_SIZE_FILTER) $(MERIPPER_WINDOW_SIZE_FILTER_OPTIONS)

%.006-windows-merged-filtered-compl.bed: %.005-windows-merged-filtered.bed
	$(COMPLEMENT_BED) $(COMPLEMENT_BED_OPTIONS) -i $^ > $@

# PER SAMPLE
$(foreach S, $(SAMPLES), $(eval $(call SAMPLE_TO_REPLICATE, $(S))))

%.peaks-compl.bed:
	$(CAT) $^ | sort -k1,1 -k2,2n | $(MERGE_BED) > $@

%.peaks.bed: %.peaks-compl.bed
	$(COMPLEMENT_BED) $(COMPLEMENT_BED_OPTIONS) -i $^ > $@

%.peaks-split.bed: %.peaks.bed
	$(MERIPPER_WINDOWSPLITTER) $(MERIPPER_WINDOWSPLITTER_OPTIONS)

# CUFFLINKS
%.cufflinks.transcripts.gtf: %-cufflinks.transcripts.gtf
	ln -sf $*-cufflinks/transcripts.gtf $*.cufflinks.transcripts.gtf

%.cufflinks.isoforms.fpkm_tracking: %-cufflinks/isoforms.fpkm_tracking
	ln -sf $*-cufflinks/isoforms.fpkm_tracking $*.cufflinks.isoforms.fpkm_tracking

%.cufflinks.genes.fpkm_tracking: %-cufflinks/genes.fpkm_tracking
	Ln -sf $*-cufflinks/genes.fpkm_tracking $*.cufflinks.genes.fpkm_tracking

%-cufflinks.transcripts.gtf %-cufflinks/isoforms.fpkm_tracking %-cufflinks/genes.fpkm_tracking: %.filt.sorted.bam
	$(CUFFLINKS) $(CUFFLINKS_OPTIONS)

# REFSEQ ANNOTATION
%.refseq.exons.bed: %.filt.bam $(REFSEQ_EXONS)
	$(COVERAGEBED) -split -abam $< -b $(REFSEQ_EXONS) | sort -k1,1 -k4,4 > $@

#%.refseq.genes.bed: %.filt.bam $(REFSEQ_GENES)
#	$(COVERAGEBED) -split -abam $< -b $(REFSEQ_GENES) | sort -k1,1 -k4,4 > $@

%.refseq.genes.rpkm.bed: $(REFSEQ_GENES) %.refseq.exons.bed %.filt.bam.rc.txt
	$(PERL) $(RPKM_PL) <(sort -k1,1 -k4,4 $<) $(wordlist 2,3,$^) > $@

# rpkm
$(ALIGNER).rpkm.txt: $(TMP_DIR)/$(GENOME).refseq.genes.cut $(addprefix $(TMP_DIR)/rpkm-, $(addsuffix .cut, $(REFSEQ_GENES_RPKM_TARGETS)))
	echo "chr txStart txEnd name score strand cdsStart cdsEnd blank exonCount exonStarts exonLengths $(FILES)" | sed 's/[ ]\+/\t/g' > $@
	paste $^ >> $@

$(TMP_DIR)/rpkm-%.refseq.genes.rpkm.bed.cut: %.refseq.genes.rpkm.bed | $(TMP_DIR)
	cut -f14 $< > $@

$(TMP_DIR)/$(GENOME).refseq.genes.cut: $(REFSEQ_GENES) | $(TMP_DIR)
	cut -f1-12 <(sort -k1,1 -k4,4 $<) > $@

# counts
$(foreach peaks, $(PEAKS_SPLIT), $(eval $(addsuffix .counts.txt, $(peaks)): $(peaks) $(addprefix $(PEAKS_COUNTS_DIR)/$(peaks)., $(addsuffix .$(ALIGNER).filt.bam.counts.cut, $(FILES)))))
$(foreach peaks, $(PEAKS_SPLIT), $(foreach file, $(FILES), $(eval $(addprefix $(PEAKS_COUNTS_DIR)/$(peaks)., $(addsuffix .$(ALIGNER).filt.bam.counts.cut, $(file))): $(peaks) $(file).$(ALIGNER).filt.bam | $(PEAKS_COUNTS_DIR) ) ) )

$(foreach peaks, $(PEAKS), $(eval $(addsuffix .counts.txt, $(peaks)): $(peaks) $(addprefix $(PEAKS_COUNTS_DIR)/$(peaks)., $(addsuffix .$(ALIGNER).filt.bam.counts.cut, $(FILES)))))
$(foreach peaks, $(PEAKS), $(foreach file, $(FILES), $(eval $(addprefix $(PEAKS_COUNTS_DIR)/$(peaks)., $(addsuffix .$(ALIGNER).filt.bam.counts.cut, $(file))): $(peaks) $(file).$(ALIGNER).filt.bam | $(PEAKS_COUNTS_DIR) ) ) )

$(ALIGNER).read-counts.txt: $(TMP_DIR)/$(GENOME).refseq.genes.cut $(addprefix $(TMP_DIR)/rc-, $(addsuffix .cut, $(REFSEQ_GENES_RPKM_TARGETS)))
	echo "chr txStart txEnd name score strand cdsStart cdsEnd blank exonCount exonStarts exonLengths $(FILES)" | sed 's/[ ]\+/\t/g' > $@
	paste $^ >> $@

$(TMP_DIR)/rc-%.refseq.genes.rpkm.bed.cut: %.refseq.genes.rpkm.bed | $(TMP_DIR)
	cut -f13 $< > $@

$(ALIGNER).mapped.txt: $(addsuffix .$(ALIGNER).filt.bam.rc.txt, $(FILES))
	echo "$(FILES)" | sed 's/[ ]\+/\t/g' > $@
	paste $^ >> $@

$(ALIGNER).library-size.txt: $(addsuffix .fastq.gz.rc.txt, $(FILES))
	echo "$(FILES)" | sed 's/[ ]\+/\t/g' > $@
	paste $^ >> $@

# ANNOTATIONS
%.bam.rc.txt: %.bam
	$(SAMTOOLS) view -c -F 0x4 $^ > $@

%.fastq.gz.rc.txt: %.fastq.gz
	expr `$(ZCAT) $^ | wc -l` / 4 > $@ 

%.$(ALIGNER).bam.frac.txt: %.$(ALIGNER).bam.rc.txt %.fastq.gz.rc.txt
	$(PERL) $(FRAC_PL) `$(CAT) $^` > $@

%.$(ALIGNER).filt.bam.frac.txt: %.$(ALIGNER).filt.bam.rc.txt %.fastq.gz.rc.txt
	$(PERL) $(FRAC_PL) `$(CAT) $^` > $@

$(TMP_DIR) $(PEAKS_COUNTS_DIR):
	mkdir -p $@

%.$(ALIGNER).peaks-split.bed.tpp.n0.pdf: %.$(ALIGNER).peaks-split.bed | $(TMP_DIR)
	$(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT) $(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT_OPTIONS) --N 0

%.$(ALIGNER).peaks-split.bed.tpp.n1.pdf: %.$(ALIGNER).peaks-split.bed | $(TMP_DIR)
	$(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT) $(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT_OPTIONS) --N 1

%.$(ALIGNER).peaks-split.bed.tpp.n2.pdf: %.$(ALIGNER).peaks-split.bed | $(TMP_DIR)
	$(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT) $(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT_OPTIONS) --N 2

$(foreach S, $(SAMPLES), $(eval $(S).$(ALIGNER).tpp.n0.pdf: $(addsuffix .$(ALIGNER).005-windows-merged-filtered.bed, $($(S))) ))
$(foreach S, $(SAMPLES), $(eval $(S).$(ALIGNER).tpp.n1.pdf: $(addsuffix .$(ALIGNER).005-windows-merged-filtered.bed, $($(S))) ))
$(foreach S, $(SAMPLES), $(eval $(S).$(ALIGNER).tpp.n2.pdf: $(addsuffix .$(ALIGNER).005-windows-merged-filtered.bed, $($(S))) ))

$(DIRNAME).$(ALIGNER).tpp.n0.pdf: $(addsuffix .$(ALIGNER).peaks-split.bed, $(SAMPLES))
	$(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT) $(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT_OPTIONS) --N 0

$(DIRNAME).$(ALIGNER).tpp.n1.pdf: $(addsuffix .$(ALIGNER).peaks-split.bed, $(SAMPLES))
	$(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT) $(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT_OPTIONS) --N 1

$(DIRNAME).$(ALIGNER).tpp.n2.pdf: $(addsuffix .$(ALIGNER).peaks-split.bed, $(SAMPLES))
	$(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT) $(MERIPPER_TRANSCRIPTOME_PROFILE_PLOT_OPTIONS) --N 2


# PEAK COUNTING
$(PEAKS_COUNTS_DIR)/%.cut:
	$(COVERAGEBED) -split -abam $(word 2, $^) -b $< -counts | sort -k1,1 -k2,2n | cut -f4 > $@

%.peaks-split.bed.counts.txt %.peaks.bed.counts.txt: 
	echo \#chr start end $(FILES) | sed 's/ /\t/g' > $@
	paste $+ >> $@

clean:
	rm *.bam *.sam *.sai *.bed *.pdf
