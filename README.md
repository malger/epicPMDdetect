# epicPMDdetect

Software that detects PMDs (Partially Methylated Domains) from Illumina (EPIC) Infinium Methylation Assays. Based on MethylseekR, adapted to support EPIC assays, using KNN datapoint selection 
and alpha-value smoothing.

### About PMDs
* PMDs are regions of intermediate methylation, that are highly disordered between CpGs.
 (in contrast the background methylation is polarized: highly methylated >=70% or depleted ) [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6161375/)
* PMDs are celltype-specific [link](https://www.biorxiv.org/content/10.1101/249334v1.full)
* PMDs can be related to closed chromatin compartments [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574526/)

### Method

The local methylation level distribution is approximated using a beta-distribution.
Its α,β parameter determines the shape. Here both parameters are set equal, meaning α=β.
Values that fullfill α > 1, result in a probability density function shape that favors 
intermediate methylation levels (PMDs). Values that fullfill α < 1, will result in a
probability density function shape that favors the polarized methylation levels (background)

Each covered CpG is assigned to an estimated α value by using its neighboring CpGs to determin the
most likely value. Those estimated alpha values are used as input to an 2-state Hidden Markov Model (HMM)
that predictes whenever a change between PMD/BG is likely. A solution (segmentation) to the HMM is
found using the Viterbi algorithm.

This method is adapted from [**MethylseekR (MSR)**](http://nar.oxfordjournals.org/content/early/2013/07/04/nar.gkt599.long):

> "Identification of active regulatory regions from DNA methylation data"
> Lukas Burger and Dimos Gaidatzis and Dirk Schubeler and Michael B. Stadler
> 2013 Nucleic Acids Research


Also [**RnBeads**](https://rnbeads.org) is used for EPIC data processing and representation 

> "Compehensive Analysis of DNA Methylation Data with RnBeads"
> "Yassen Assenov and Fabian Mueller and Pavlo Lutsik and Joern Walter and Thomas Lengauer and Christoph Bock"
> 2014 Nature Methods 11

### Adjustments to work on EPIC data

#### Support EPIC data

MSR was developed with WGBS (Whole genome bisulfite sequencing) data in mind. 
Therefore it operates on read counts, comparing reads of methylated CpGs to the total number of reads.
EPIC assays use [red/green light intensities](https://en.wikipedia.org/wiki/Illumina_Methylation_Assay) to determin the methylation level. The total number of reads is calculated by the sum of light intensities, light intensities covering methylated CpGs are passed as methylated reads.

#### Respect difference in data point distributions

In comparison to WGBS, EPIC assays only cover a small subset of CpGs from the human genome. Therefore 
size constrains (distance cutoff) and nearest neighbor selection using KNN ensure that only data points in reasonable distance are used for alpha value estimation. 

#### Reduce influence of certain parameter setups

Per default different values for k (neighboorhood size) and distance cutoff are evaluated to estimate alpha. Those alpha value estimations are averaged to diminish the effect of individual parameter setups 
to alpha.

#### Mark unreliable segments 

Some areas have too few data points to reliably assign segments. Those are marked accordingly

### How to use

load the library

```R
library(epicPMDdetect)
```

generate a csv sampleSheet that describes your data, herefore refer to the instructions in [RnBeads Manual](https://bioconductor.org/packages/release/bioc/vignettes/RnBeads/inst/doc/RnBeads.pdf).

Now load the data using rnbeads by calling *readEPIC_idats*.
The *sample.col.name* musst be the column name of the csv holding the sample names 


```R
data = readEPIC_idats(sample.sheet,idat.dir,sample.col.name,preprocess=T)
```

Specific which samples should be processed (default all) by providing an vector of names 
(optional a settings object can be passed that controls further parameters.)

```R
segmentRnbSet(data,outputFolder,samples=NULL)

```

The outputFolder contains the segmentations for each specified sample.
They are named:

**[sampleName].seg.bed.gz**




