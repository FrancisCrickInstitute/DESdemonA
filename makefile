tmp1 =  $(lastword $(sort $(wildcard R-*-local))) # highest version number.
init_R = $(strip $(tmp1)) # remove whitespace
labscipro = $(subst /, ,$(subst /camp/stp/babs/working/kellyg/projects,,$(CURDIR)))
lab = $(word 1,$(labscipro))
sci = $(word 2,$(labscipro))
firstname = $(word 1, $(subst ., ,$(sci)))
pro = $(word 3,$(labscipro))

.PHONY:  clean dirs links init rnaseq analysis
################################################################
# Generic Project recipes
################################################################
init: links dirs
clean:
	@rm -f analyse.r.out.log
	@rm -f analyse.r.err.log
	@rm -f nohup.out

links:
	ln -sf /camp/stp/babs/www/kellyg/public_html/LIVE/projects/$(lab)/$(sci)/$(pro)/ www
	ln -sf /camp/stp/babs/outputs/$(lab)/$(sci)/gavin.kelly/$(pro)/ outputs

dirs: 
	mkdir -p results
	mkdir -p objects
	mkdir -p data

rnaseq:
	module load nextflow/0.30.0 ;\
	nohup nextflow run -w work rnaseq.nf -params-file align.yml &

analysis: analyse.r
ifeq ($(executor),slurm)
	[[ -z "$(init_R)" ]] || source "$(init_R)" ;\
	sbatch --wrap="Rscript $<" -J "$(firstname)_$<" -e "$<.err.log"  -o "$<.out.log"
else
	[[ -z "$(init_R)" ]] || source "$(init_R)" ;\
        Rscript $<
endif

