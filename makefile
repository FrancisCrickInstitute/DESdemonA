me = gavin.kelly
pathwords = $(subst /, ,$(CURDIR)))
sci = $(word 8,$(pathwords))
www = $(subst /babs/working/${USER}/,/babs/www/${USER}/public_html/LIVE/,$(CURDIR))
scratch = $(subst /babs/working/,/babs/scratch/,$(CURDIR))
outputs = $(subst /$(sci)/,/$(sci)/$(me)/,$(subst working/${USER}/projects,outputs,$(CURDIR)))


.PHONY:  clean dirs links init rnaseq analysis
################################################################
# Generic Project recipes
################################################################
init: links dirs
clean:
	@rm -f analyse.r.out.log
	@rm -f analyse.r.err.log
	@rm -f nohup.out

init: dirs links

links:
	mkdir -p $(www); ln -sf $(www) www
	mkdir -p $(outputs); ln -sf $(outputs) outputs
	mkdir -p $(scratch); ln -sf $(scratch) scratch

dirs: 
	mkdir -p results
	mkdir -p objects
	mkdir -p data

rnaseq:
	module load nextflow/0.30.0 ;\
	nohup nextflow run -w scratch/work -resume rnaseq.nf -params-file align.yml &

analysis: init_R =  $(lastword $(sort $(wildcard R-*-local)))# highest version number.
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


stabilise_links : dir=results
stabilise_links : iname=*
stabilise_links:
	for f in $$(find $(dir) -iname "$(iname)" -type l); do \
	mv "$$(readlink -e $$f)" "$$f";\
	done
