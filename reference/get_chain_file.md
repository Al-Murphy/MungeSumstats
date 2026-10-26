# Download chain file for liftover

Download chain file for liftover

## Usage

``` r
get_chain_file(
  from = c("hg38", "hg19"),
  to = c("hg19", "hg38"),
  chain_source = c("ucsc", "ensembl"),
  save_dir = tempdir(),
  verbose = TRUE
)
```

## Source

[UCSC chain
files](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)

[Ensembl chain
files](https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/)

## Arguments

- from:

  genome build converted from ("hg38", "hg19")

- to:

  genome build converted to ("hg19", "hg38")

- chain_source:

  chain file source used ("ucsc" as default, or "ensembl")

- save_dir:

  where is the chain file saved? Default is a temp directory

- verbose:

  extra messages printed? Default is TRUE

## Value

loaded chain file for liftover
