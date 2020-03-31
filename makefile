################################################################
#### RNASeq recipes
################################################################

.PHONY:   rnaseq analysis

rnaseq:
	module load nextflow/0.30.2 ;\
	nohup nextflow run -w scratch/work \
		-params-file params.yml \
		-r v0.1.0 \
		--results_dir data \
		http://github.com/crickbabs/BABS-RNASeq  &

analysis: init_R = $(lastword $(sort $(wildcard R-*-local)))# highest version number.
analysis: analyse.r
	[[ -z "$(init_R)" ]] || source "$(init_R)" ;\
ifeq ($(executor),slurm)
	sbatch --wrap="Rscript $<" -J "$(sci)_$<" -e "$<.err.log"  -o "$<.out.log"
else
        Rscript $<
endif


fastq_folder = /camp/stp/sequencing/inputs/instruments/data/$(lab)/$(sci)/$(asf)/primary_data/
paired = $(shell find -L ${fastq_folder} -path "*/fastq/*" -name "*_R2_*.fastq.gz" | wc -l)
extdata/design.csv:
ifeq (${paired}, 0)
	find -L ${fastq_folder} -path "*/fastq/*" -name "*_R1_*.fastq.gz" -printf "%P\n" | awk 'BEGIN { FS = "[_/]" ; OFS="," ; print "file,date,machine,run,sample,id,lane" } { print "${fastq_folder}"$$0, $$1, $$2, $$3, $$6, $$7, $$8 }' > $@
else
	find -L ${fastq_folder} -path "*/fastq/*" -name "*_R1_*.fastq.gz" -printf "%P\n" | awk 'BEGIN { FS = "[_/]" ; OFS="," ; print "file1,file2,date,machine,run,sample,id,lane" } {r2=$$0; gsub("_R1_","_R2_",r2); print "${fastq_folder}"$$0, ${fastq_folder} r2,$$1, $$2, $$3, $$6, $$7, $$8}' >@
endif
