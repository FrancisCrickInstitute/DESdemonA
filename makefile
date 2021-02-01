R := R-4.0.2-BABS
res_dir := results


ifeq (,$(wildcard .git))
TAG := unversioned
VERSION:= unversioned
else
TAG := $(shell git describe --tags --dirty --always --long)
VERSION := $(shell git describe --tags --abbrev=0)
endif

analysis : contained := FALSE
analysis: 01_analyse.r  data/rsem_dds.rda version_dir R.bib
	$(R) -e "rmarkdown::render('01_analyse.r',\
	  output_file='output_$(TAG).html',\
	  output_options = list(self_contained=$(contained)),\
	  params=list(res_dir='$(res_dir)/$(VERSION)'))"
	rm -f $(res_dir)/$(VERSION)/output_$(TAG).html
	mv output_$(TAG).html $(res_dir)/$(VERSION)
	rm -rf $(res_dir)/$(VERSION)/output_$(TAG)_files
	mv output_$(TAG)_files $(res_dir)/$(VERSION)


data/rsem_dds.rda: init.r $(wildcard inst/extdata/rsem/*.genes.results)
	$(R) -e "source('init.r')"

version_dir:
	mkdir -p $(res_dir)/$(VERSION)

R.bib: analyse.r
	$(R) -e "pd <- getParseData(parse('$<', keep.source=TRUE));\
	libreq <- pd\$$text[pd\$$line1 %in% pd\$$line1[pd\$$text=='library' | pd\$$text=='require'] & pd\$$token=='SYMBOL'];\
	libreq <- unique(c('base', libreq, pd\$$text[pd\$$token=='SYMBOL_PACKAGE']));\
	knitr::write_bib(libreq, file='$@')"

design.csv:
	ls asf/fastq/*_R1_*fastq.gz |
	awk  'BEGIN {FS = "[_/]"; print "sample,file1,file2"} ; {r2=$$0;  sub(/_R1_/, "_R2_", r2); print $$3","$$0","r2}' > design.csv

alignment: design_nf.csv
	ml purge ;\
	ml Nextflow/20.12.0-edge ;\
	ml Singularity/3.4.2 ;\
	ml CAMP_proxy ;\
	export NXF_SINGULARITY_CACHEDIR=/camp/apps/misc/stp/babs/nf-core/singularity/rnaseq/3.0/ ;\
	nextflow run nf-core/rnaseq \
	--input $< \
	--outdir results \
	--aligner star_rsem \
	--fasta /camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome/Mus_musculus.GRCm38.dna_sm.toplevel.fa \
	--gtf /camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf \
	--rsem_index /camp/stp/babs/working/patelh/genome/GRCm38/release-95/index/rsem/star-2.7.6a/ \
	-w $(shell readlink -f scratch)/nf_work \
	-profile crick \
	-c custom.config \
	-resume \
	-with-tower \
	-r 3.0

-include babs.mk
