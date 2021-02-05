## DESDemonA Philosophy

This is a wrapper round DESeq2 which standardises the process of
generating `DESeq2::results` objects for a specifiable set of
comparisons. Each comparison is done in the context of a model (a
specification of what covariates need to be accounted for in
predicting the expression of any transcript). And each model is
formulated in a sample-set (which is usually all the samples, but can
be changed to drop samples, or look at a subset in an isolated way).

There are two inputs to the analysis pipeline: a raw DESeq2 object;
and an analysis specification.  The DESeq2 object should be stored as
`data/rsem_dds.rda` which contains all the counts in its `counts`
assay, and all the predictors should be stored as columns of the
`colData` slot. We'll show a simple way of turning the results of the
NF-Core RNASeq pipeline into a suitable object later.

All the hints as to what analyses are to be carried out are stored in
a `.spec` file, and there can be as many of these as you want - so you
could have an `overview.spec` file that tried out a handful of
approaches, and a `publication.spec` that recorded the chosen
approach. Each `.spec` file is in fact an R script (and so should be
treated with as much caution, security-wise, as any other script), but
is expected to be of a certain form, the first criterion is that it
should start and end as follows:

```
specification(
...
)
```

### Sample sets
The highest level of any DESDemonA analysis is a choice as to what
samples are to be analysed together. Often this is trivial, and the
only sensible analysis is to look at all the samples together. But
even in the simplest experiment, it might be observed for example that
one sample is an outlier and we'd like to examine the effect of
omitting or retaining this sample. We can accomplish this easily:

```
specicifation(
  sample_sets =list( 
    all = sample_set(
      subset = TRUE
    ),
	no_outlier=sample_set(
	  subset = sample_id != "id_of_bad_sample"
    )
  )
)
```

So a collection of samples is called a `sample_set`, and we have two
(a `list`) of theme here. The `subset` expression of a sample\_set is
evaluated in the context of the `colData`, so the first one always
evaluates to `TRUE`, denoted that we include all samples; the second
one evaluates to `FALSE` for a particular sample (assuming we have a
sample\_id column in our `colData`), dropping that sample.

Other reasons we may want to subset our samples into different
configurations is: if there are two very different cohorts of samples
where transcripts are likely to have substantially different noise
profiles and we want to avoid contaminating one set of samples with
the estimates of dispersion that hold for the other set. Of course it
then becomes impossible to _statistically test_ between the two
cohorts, but we can have a third set that still combines them:

```
specicifation(
  sample_sets =list( 
    all = sample_set(
      subset = TRUE
    ),
	tumour=sample_set(
	  subset = condition == "tumour"
    ),
	normal=sample_set(
	  subset = condition == "normal"
    )

  )
)
```

It is possible to transform some of the predictors (columns of
`colData`) at this stage, but we will (describe this)[later]. But for
now we'll look at the main analysis attribute of a sample set, namely
the model we'll use to estimate the expression in the various
experimental conditions those samples were taken from:

### Models

DESeq2 uses R's standard one-sided formulas to inform the analysis as
to which predictors to use in fitting a negative-binomial model of
expression counts. The `colData` slot should contain all the potential
predictors, and the `design` slot of a DESeq2 object encodes a
considered choice of which are relevant to the biological
question. Obviously, if the question is to find which genes are
differentially expressed between two treatments, then a predictor
annotating which treatment a sample is subject to is necessary. But
quite often there will be other covariates need to be accounted for,
such as batch; or perhaps there's a subtler question of finding genes
that respond differently to treatment depending on the genotype of the
sample.

We may want to consider alternative models, for instance whether or
not to include a batch effect: if the batch effect is negligible we
won't want to use up degrees of freedom (ie power) estimating it, but
if it's substantial then we want to make our estimates more precise by
accounting for it. So we allow, in any given sample_set, for there to
be multiple models:

```
specicifation(
  sample_sets =list( 
    all = sample_set(
      subset = TRUE,
	  models = list(
	    accurate = model(
		  design = ~treatment + batch
		),
	    powerful = model(
		  design = ~treatment
		)
	  )
    )
  )
)
```

So once again the plural 'models' is a list of individual `model`
objects, which for the moment we've just specified by their `design`,
exactly as we would for a `DESeq` object. Once we have a model, we are
then in a position to test the significance of contrasts (e.g.
individual log fold-changes between condtions) and terms in the model
(e.g whether batch effect is necessary or not):

### Comparisons

DESeq2 has multiple ways of specifying a contrast. All of these are
supported, and again as 'comparisons' is plural, we can have a list of
them, and have DESDemonA loop through all of them for the parent
model:

```
specicifation(
  sample_sets =list( 
    all = sample_set(
      subset = TRUE,
	  models = list(
	    accurate = model(
		  design = ~treatment + batch,
		  models = list(
			comp1 = "treat_vs_ctrl",
			comp2 = c("treatment", "treat", "ctrl"),
			comp3 = list(c("treat_vs_ctrl"), c("vehicle_vs_ctrl"), listValues=c(1,-1))
			comp4 = c(1,0,-0.5,-0.5)
		  )
		)
	  )
    )
  )
)
```

So all of the traditional ways of specifying a specific contrst are
represented, respecively:

- by a single character, traditionally done in DESeq2 with
  `results(..., name="treat_vs_ctrl")`
  
- A character triplet

- A two-component list of characters indicating the numerator and
  denominator (with possible `listValues` set)
  
- A numeric vector.

You are referred to the DESeq2 reference manual for further details,
but one advantage of DESDemonA is that it allows a simple enumeration
of complex comparisons using machinery from the `emmeans` package:

#### Automatic comparison enumeration

We often find ourselves copying a contrast multiple times, with slight
changes to elicit comparisons between different conditions. For
example if our samples have been subject to three different treatment
regimes, and a control, so treatment={_control_, _vehicle_,
_standard_, _novel\_treatment_}, and we wanted to compare everything to
_control_, then we'd need to hand-write three comparisons. Instead we
can now write

```
specicifation(
  sample_sets =list( 
    all = sample_set(
      subset = TRUE,
	  models = list(
	    accurate = model(
		  design = ~treatment + batch,
		  models = list(
			mult_comp(trt.vs.ctrl ~ treatment, ref="control")
		  )
		)
	  )
    )
  )
)
```

where the `trt.vs.ctrl` is a 'keyword' instructing as to which comparisons to
carry out - in this case every value of `treatment` against the ones
labelled "control". There are a number of other helpful keywords:

- `rev.pairwise` is the most 'verbose', in that it will generate every
  pair of comparisons, in this case there are six combinations, and it
  can go up rapidly. Due to a slightly annoying convention, `pairwise`
  choses an unintuitive ordering that would result in comparisons such
  as 'control vs novel\_treatment' (and a positive fold change would mean
  higher expression in the control than the novel treatment), and
  `rev.pairwise` has the more conventional ordering of
  'novel\_treatment vs control' and it's more intuitive positive log
  fold-change indicating higher expression in the treated to the
  control. This assumes the levels of the factor are given in the
  conventional order of the previous paragraph.
  
- `consec` takes adjacent pairs along the levels, so again quite a
  natural 'vehicle vs control', 'standard vs vehicle' and finally
  'novel\_treatment vs standard'
  
- `trt.vs.ctrl` choses one baseline (which we can select with an
  additional `ref` argument, rather than having to `relevel` the data)
  and compares everything else against it.
  
When we don't have an interaction term in our model, this is
everything we'd ever need, as for example the difference between
'novel\_treatment' and 'control' is by construction independent of
which batch we're in.

#### Interaction terms

When our model includes an interaction between two main effects, say
between treatment and genotype, then the magnitude of difference
between two treatments can depend on the genotype. This leads to two
types of follow-up question: stratifying the data by one factor, and
then examining the effect of the other factor in each individual
stratum, or; examining the interaction itself, so looking to see if
the magnitude of effect is different in one stratume than another.

To make this more concrete for the genotype x treatment case, we could
ask: which genes are responding to treatment in the KO; which genes
are responding to treatment in the WT; which genes are differential
between WT and KO in the untreated. These are all questions of the
first type, and we can automatically generate their contrasts by:

```
specicifation(
  sample_sets =list( 
    all = sample_set(
      subset = TRUE,
	  models = list(
	    accurate = model(
		  design = ~treatment *  genotype,
		  models = list(
			mult_comp(rev.pairwise ~ treatment | genotype)
			mult_comp(rev.pairwise ~ genotype  | treatment)
		  )
		)
	  )
    )
  )
)
```

so the `|` symbol is to be read as "stratified by" (or "conditioned
on"). The first `mult_comp` is stratifying the treatment, and will
generate the pairwise comparisons within treatment, separately for
each different genotype. The second is the dual quesiton: For each
treatment, find the genes that are 'markers of genotype' within that
treatment group.

The second type of question, where we want to investigate a
'difference of differences', is achieved by the following grammar:

```
mult_comp(rev.pairwise+rev.pairwise ~ treatment + genotype, interaction=TRUE)
```

The terms on the right hand side of the formula enumerate the
predictors we're concerned in (there's no need to have their
interaction - that's already known because of the `design`), and the
left hand side specifies the contrasts that are going to be considered
(so any of the keywords from the previous section) - there need's to
be either one (which will be applied to all predictors), or one per
factor on the other side of the formula (as here) in the respective
order.

The above is probably quite intimidating, so please ask for help on
this - I'll try to write something a bit more friendly.


### Rank deficiency 

One of the most common hurdles in complex designs s the
rank-deficiency problem. The way the experiment has been described can
result in the analysis struggling to identify an estimate for a
particular experimental condition. This happens for one of three
reasons:

- The design has a fatal flaw. If all your treated samples are in one
  batch, and all your control samples in another batch, you cannot
  estimate the extent on how these individually influence expression.
  
- Some conditions have no samples. For example, in our genotype x
  treatment experiment, we might not have a 'vehicle' treatment in one
  of our genotypes.
  
- There is a nesting indicative of repeated measurements. Individuals
  might have submitted multiple samples (perhaps a time-course), and
  also be members of larger groups (their disease status, which is
  assumed constant throughout the time-course).
  
The latter two are remediable through some tweaks to the inputs of
DESeq2, the first will need the advice of a statistician to see if
anything is recoverable and at what compromise.

The 'missing condition' case is easily remedied by supplying the
option `drop\_unsupported\_combinations=TRUE` to any _model_ you want
this to apply to. The 'nesting' case normally requires a tweak to the
way factors are coded, and we have made it easy to achieve this
through a `transform` option to any _sample\_set_.  Both options are
illustrated below:

```
specicifation(
  sample_sets =list( 
    all = sample_set(
      subset = TRUE,
	  transform = mutate(unit = recode_within(unit, genotype))
	  models = list(
	    accurate = model(
		  design = ~timepoint *  genotype + person:genotype,
		  models = list(
			mult_comp(rev.pairwise ~ timepoint | genotype),
		  drop_unsupported_combinations=TRUE
		  )
		)
	  )
    )
  )
)
```

The 'missing condition' case is achieved by just the one line; the
nesting is more complex - we firstly have to have a 'unit' column in
our `colData` that identifies which biological unit (person, animal,
plate) a sample belongs to. Then the `transform` line calls a `recode\_
within` function that tells us that unit is nested within genotype
here. And we alter the `design` in accordance with the DESeq2 vignette
on within- and between- group comparisons.

The `transform` statement can be used to adjust any of the `colData`
columns in-line, though it might be better to generate then in the
original DESeq2 data object that is taken as input.
