TAG := $(shell git describe --tags --dirty --always --long)
VERSION := $(shell git describe --tags --abbrev=0)
R := R-3.6.0-local

analysis : contained := FALSE
analysis: analyse.r  data/rsem_dds.rda version_dir
	$(R) -e "rmarkdown::render('analyse.r',\
	  output_file='output_$(TAG).html',\
	  output_options = list(self_contained=$(contained)),\
	  params=list(res_dir='results/$(VERSION)'))"
	mv output_$(TAG).html results/$(VERSION)/
	[ ! -d output_$(TAG)_files ]  || mv output_$(TAG)_files results/$(VERSION)/

data/rsem_dds.rda: init.r $(wildcard inst/extdata/rsem/*.genes.results)
	$(R) -e "source('init.r')"

version_dir:
	mkdir -p results/$(VERSION)

R.bib: analyse.r
	$(R) -e "pd <- getParseData(parse('analyse.r', keep.source=TRUE));\
	libreq <- pd\$$text[pd\$$line1 %in% pd\$$line1[pd\$$text=='library' | pd\$$text=='require'] & pd\$$token=='SYMBOL'];\
	libreq <- unique(c('base', libreq, pd\$$text[pd\$$token=='SYMBOL_PACKAGE']));\
	knitr::write_bib(libreq, file='$@')"

design.csv:
	ls asf/fastq/*_R1_*fastq.gz |
	awk  'BEGIN {FS = "[_/]"; print "sample,file1,file2"} ; {r2=$$0;  sub(/_R1_/, "_R2_", r2); print $$3","$$0","r2}' > design.csv

alignment: params.yml
	module load nextflow/0.30.2 ;\
	nohup nextflow run -w scratch/work \
	-params-file params.yml \
	-resume \
	-latest \
	-r 9e35296 \
	--results_dir data \
	http://github.com/crickbabs/BABS-RNASeq  &


#### Remainder is to initialise projects in a Crick-specific way
#### There should be no need to run this again
# EMAIL environment variable should be set, so we can work out correct 'outputs' directory.
me := $(firstword $(subst @, ,${EMAIL}))
babsfield = $(shell awk '$$1 == "$(1):"{for (i=2; i<=NF; i++) print $$i}' .babs)
lab := $(call babsfield,Lab)
sci := $(shell echo $(firstword $(subst @, ,$(call babsfield,Scientist))) | tr A-Z a-z)
type := $(call babsfield,Type)
lims := $(call babsfield,Lims)
# albert.einstein => Albert Einstein
scientist := $(shell echo $(sci) | sed --expression="s/\b\(.\)/\u\1/g; s/\./ /")
# for R package author albert.einstein => "Albert","Einstein note no final quote - it's added explicitly
person := $(shell echo $(me) | sed --expression="s/\b\(.\)/\"\u\1/g; s/\./,/")
Me := $(shell echo $(me) | sed --expression="s/\b\(.\)/\u\1/g; s/\./ /")
strproject := $(call babsfield,Project)
project := $(shell echo $(strproject) | sed --expression="s/[^a-zA-Z0-9]/_/g")
www := /camp/stp/babs/www/${USER}/public_html/LIVE/projects/$(lab)/$(sci)/$(project)
scratch := /camp/stp/babs/scratch/${USER}/$(lab)/$(sci)/$(project)
outputs := /camp/stp/babs/outputs/$(lab)/$(sci)/$(me)/$(project)

config:
	mkdir -p data objects results inst/extdata
	mkdir -p $(www) $(scratch) $(outputs)
	ln -sfn $(www) www
	ln -sfn $(scratch) scratch
	ln -sfn $(outputs) outputs
ifneq ($(lims),{{lims}})
	ln -sfn /camp/stp/sequencing/inputs/instruments/data/$(lab)/$(sci)/$(lims)/primary_data/ asf
	cp -n /camp/stp/sequencing/inputs/instruments/data/$(lab)/$(sci)/$(lims)/$(lims)_design.csv inst/extdata/design.csv
endif
	sed -i 's/{{descrip}}/Analysis for $(scientist) in $(lab) lab/g' DESCRIPTION
	sed -i 's/{{version}}/$(VERSION)/g' DESCRIPTION
	sed -i 's/{{email}}/${EMAIL}/g' DESCRIPTION
	sed -i 's/{{person}}/${person}"/g' DESCRIPTION
	for r in *.{r,R,rmd,Rmd} DESCRIPTION ; do \
	sed -i 's/{{project}}/$(strproject)/g' $$r ; \
	sed -i 's/{{package}}/babs$(type)/g' $$r ; \
	sed -i 's/{{author}}/$(Me)/g' $$r ; \
	done
	git commit -a -m "Standard starting point"
	git tag v0.0.1

