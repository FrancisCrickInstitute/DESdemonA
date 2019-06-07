include ../generic/makefile

.PHONY:   rnaseq analysis

rnaseq:
	module load nextflow/0.30.0 ;\
	nohup nextflow run -w scratch/work -resume rnaseq.nf -params-file align.yml &

analysis: init_R = $(lastword $(sort $(wildcard R-*-local)))# highest version number.
analysis: analyse.r
	[[ -z "$(init_R)" ]] || source "$(init_R)" ;\
ifeq ($(executor),slurm)
	sbatch --wrap="Rscript $<" -J "$(sci)_$<" -e "$<.err.log"  -o "$<.out.log"
else
        Rscript $<
endif

data/design.csv: paired = $(shell ls fastq/*_R2_* | wc -l)
data/design.csv: fastq
ifeq ($paired, 0)
	ls fastq | grep "_R1_" | awk 'BEGIN { FS = "_" ; OFS="," ; print "file,sample,id,lane" } { print "fastq/"$$0, $$1, $$2, $$3 }' > data/design.csv
else
	ls fastq | grep "_R1_" | awk 'BEGIN { FS = "_" ; OFS="," ; print "file1,file2,sample,id,lane" } { r2=$$0; gsub("_R1_","_R2_",r2); print "fastq/"$$0,"fastq/" r2 , $$1, $$2, $$3 }' > data/design.csv
endif

