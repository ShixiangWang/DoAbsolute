
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DoAbsolute

The goal of DoAbsolute is to automate ABSOLUTE calling for multiple
samples in parallel way.

[ABSOLUTE](https://www.nature.com/articles/nbt.2203) is a famous
software developed by Broad Institute, however the `RunAbsolute`
function is designed for computing one sample each time and set no
default values. **DoAbsolute** help user set default parameters
according to [ABSOLUTE
documentation](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/ABSOLUTE),
provide an uniform interface to input data easily and run RunAbsolute
parallelly.

More detail about how to analyze ABSOLUTE results please see [this
link](http://software.broadinstitute.org/cancer/software/genepattern/analyzing-absolute-data).

## Installation

You can install the released version of DoAbsolute with:

``` r
devtools::install_github("ShixiangWang/DoAbsolute")
```

Install ABSOLUTE, the version provided by DoAbsolute is 1.0.6. You can
find available versions at
<https://software.broadinstitute.org/cancer/cga/absolute_download>.
Users of DoAbsolute all should accept LICENCE from Firehose.

``` r
install.packages("numDeriv")
path_to_file = system.file("extdata", "ABSOLUTE_1.0.6.tar.gz", package = "DoAbsolute", mustWork = T)
install.packages(path_to_file, repos = NULL, type="source")
```

## Example

This is a basic example which shows you how to run DoAbsolute using
example data from [ABSOLUTE
documentation](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/ABSOLUTE).

Load
package.

``` r
library(DoAbsolute)
```

``` r
example_path = system.file("extdata", package = "DoAbsolute", mustWork = T)


library(data.table)
# Load Test Data ----------------------------------------------------------

# segmentation file
seg_normal =  file.path(example_path, "SNP6_blood_normal.seg.txt")
seg_solid  =  file.path(example_path, "SNP6_solid_tumor.seg.txt")
seg_metastatic  = file.path(example_path, "SNP6_metastatic_tumor.seg.txt")
# MAF file
maf_solid  = file.path(example_path, "solid_tumor.maf.txt")
maf_metastatic  = file.path(example_path, "metastatic_tumor.maf.txt")

# read data
seg_normal = fread(seg_normal)
seg_solid = fread(seg_solid)
seg_metastatic = fread(seg_metastatic)
maf_solid = fread(maf_solid)
maf_metastatic = fread(maf_metastatic)

# merge data
Seg = Reduce(rbind, list(seg_normal, seg_solid, seg_metastatic))
Maf = Reduce(rbind, list(maf_solid, maf_metastatic))

Seg$Sample = substr(Seg$Sample, 1, 15)
Maf$Tumor_Sample_Barcode = substr(Maf$Tumor_Sample_Barcode, 1, 15)

# test function
DoAbsolute(Seg = Seg, Maf = Maf, platform = "SNP_6.0", copy.num.type = "total",
           results.dir = "test", nThread = 2, keepAllResult = TRUE, verbose = TRUE)
```

## Citation

*Wang, Shixiang, et al. “APOBEC3B and APOBEC mutational signature as
potential predictive markers for immunotherapy response in non-small
cell lung cancer.” Oncogene (2018).*
