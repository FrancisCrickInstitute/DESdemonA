ifeq (,$(wildcard .babs))
babsfield = $(1)
else
babsfield = $(shell awk '$$1 == "$(1):"{for (i=2; i<=NF; i++) print $$i}' .babs)
endif

# git should be config'ed with user.email and user.name
email := $(shell git config --get user.email)
me := $(shell git config --get user.name)
lab := $(call babsfield,Lab)
sci := $(shell echo $(firstword $(subst @, ,$(call babsfield,Scientist))) | tr A-Z a-z)
type := $(call babsfield,Type)
lims := $(call babsfield,Lims)
# albert.einstein => Albert Einstein
scientist := $(shell echo $(sci) | sed --expression="s/\b\(.\)/\u\1/g; s/\./ /")
empty :=
space := $(empty) $(empty)
comma :=,
person := $(subst $(space),"$(comma)",$(me))
strproject := $(call babsfield,Project)
project := $(shell echo $(strproject) | sed --expression="s/[^a-zA-Z0-9]/_/g")
www := /camp/stp/babs/www/${USER}/public_html/LIVE/projects/$(lab)/$(sci)/$(project)
scratch := /camp/stp/babs/scratch/${USER}/$(lab)/$(sci)/$(project)
outputs := /camp/stp/babs/outputs/$(lab)/$(sci)/$(firstword $(subst @, ,${email}))/$(project)

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
	sed -i 's/{{email}}/${email}/g' DESCRIPTION
	sed -i 's/{{person}}/"${person}"/g' DESCRIPTION
	for r in *.{r,R,rmd,Rmd} DESCRIPTION ; do \
	sed -i 's/{{project}}/$(strproject)/g' $$r ; \
	sed -i 's/{{package}}/babs$(type)/g' $$r ; \
	sed -i 's/{{author}}/$(me)/g' $$r ; \
	done
	git commit -a -m "Standard starting point"
	git tag v0.0.1

