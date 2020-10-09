R := R-3.6.0-local
res_dir := results
nf_w := scratch/work
nf_params := params.yml
nf_switches := -resume -latest -r 9e35296
nf_results_dir := data
nf_url := http://github.com/crickbabs/BABS-RNASeq


ifeq (,$(wildcard .git))
TAG := unversioned
VERSION:= unversioned
else
TAG := $(shell git describe --tags --dirty --always --long)
VERSION := $(shell git describe --tags --abbrev=0)
endif

analysis : contained := FALSE
analysis: analyse.r  data/rsem_dds.rda version_dir
	$(R) -e "rmarkdown::render('analyse.r',\
	  output_file='output_$(TAG).html',\
	  output_options = list(self_contained=$(contained)),\
	  params=list(res_dir='$(res_dir)/$(VERSION)'))"
	mv output_$(TAG).html $(res_dir)/$(VERSION)/
	[ ! -d output_$(TAG)_files ]  || mv output_$(TAG)_files $(res_dir)/$(VERSION)/

data/rsem_dds.rda: init.r $(wildcard inst/extdata/rsem/*.genes.results)
	$(R) -e "source('init.r')"

version_dir:
	mkdir -p $(res_dir)/$(VERSION)

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
	nohup nextflow run -w $(nf_w) \
	-params-file $(nf_params) \
	$(nf_switches) \
	--results_dir $(nf_results_dir) \
	$(nf_url)  &

-include babs.mk
