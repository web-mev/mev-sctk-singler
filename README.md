# mev-sctk-singler

This repository contains a WDL-format Cromwell-compatible workflow for executing an auto cell typing analysis on single-cell RNA-seq data using the SingleR tool (https://www.nature.com/articles/s41590-018-0276-y) as provided through the Single-Cell Toolkit (https://github.com/compbiomed/singleCellTK).

To use, simply fill in the the `inputs.json` with the various inputs and submit to a Cromwell runner. 

Alternatively, you can pull the docker image (), start the container, and run: 

```
Rscript /opt/software/singler.R \
    -f <path to raw counts tab-delimited file> \
    -o <prefix for output file (string)> \
    --level <A (string) that specifies the granularity of cell typing. Choose from main or fine.> \
    --reference <A (string) that specifies a reference provided by SingleR. Choose from hpca, bpe, mp, dice, immgen, mouse, zeisel.> \
    --featureType <A (string) for whether to use gene symbols or Ensembl IDs when using a SingleR built-in reference. Choose from symbol, ensembl.>
```