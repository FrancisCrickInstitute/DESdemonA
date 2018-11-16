init_R =  $(lastword $(sort $(wildcard R-*-local)))# highest version number.
me = gavin.kelly
pathwords = $(subst /, ,$(CURDIR)))
sci = $(word 8,$(pathwords))
www = $(subst /babs/working/${USER}/,/babs/www/${USER}/public_html/LIVE/,$(CURDIR))
scratch = $(subst /babs/working/,/babs/scratch/,$(CURDIR))
outputs = $(subst /$(sci)/,/$(sci)/$(me)/,$(subst working/${USER}/projects,outputs,$(CURDIR)))

.PHONY:  clean init rnaseq analysis

################################################################
# Generic Project recipes
################################################################
init: dirs links

links:
	mkdir -p $(www); ln -sf $(www) www
	mkdir -p $(outputs); ln -sf $(outputs) outputs
	mkdir -p $(scratch); ln -sf $(scratch) scratch

dirs: 
	mkdir -p results
	mkdir -p objects
	mkdir -p data


clean:
	@rm -f analyse.r.out.log
	@rm -f analyse.r.err.log
	@rm -f nohup.out


################################################################
#### Project-type recipes
################################################################
rnaseq:
	module load nextflow/0.30.0 ;\
	nohup nextflow run -w work rnaseq.nf -params-file align.yml &

analysis: analyse.r
ifeq ($(executor),slurm)
	[[ -z "$(init_R)" ]] || source "$(init_R)" ;\
	sbatch --wrap="Rscript $<" -J "$(sci)_$<" -e "$<.err.log"  -o "$<.out.log"
else
	[[ -z "$(init_R)" ]] || source "$(init_R)" ;\
        Rscript $<
endif
