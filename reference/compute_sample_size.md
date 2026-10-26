# Compute (effective) sample size

Computes sample sum (as new column "N") or effective sample size (ESS)
(as new column "Neff"). Computing ESS is important as it takes into
account the proportion of cases to controls (i.e. class imbalance) so as
not to overestimate your statistical power.

## Usage

``` r
compute_sample_size(
  sumstats_dt,
  method = c("ldsc", "giant", "metal", "sum"),
  force_new = FALSE,
  append_method_name = FALSE
)
```

## Arguments

- sumstats_dt:

  Summary statistics data.table.

- method:

  Method for computing (effective) sample size.

  - "ldsc" :  
    \\Neff = (N_CAS+N_CON) \* (N_CAS/(N_CAS+N_CON)) /
    mean((N_CAS/(N_CAS+N_CON))\[(N_CAS+N_CON)==max(N_CAS+N_CON)\]))\\  
    [bulik/ldsc GitHub Issue](https://github.com/bulik/ldsc/issues/95)
    [bulik/ldsc GitHub
    code](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L321)

  - "giant" :  
    \\Neff = 2 / (1/N_CAS + 1/N_CON)\\  
    [Winkler et al. 2014,
    Nature](https://www.nature.com/articles/nprot.2014.071)

  - "metal" :  
    \\Neff = 4 / (1/N_CAS + 1/N_CON)\\  
    [Willer et al. 2010,
    Bioinformatics](https://pubmed.ncbi.nlm.nih.gov/20616382/)

  - "sum" :  
    \\N = N_CAS + N_CON\\  
    Simple summation of cases and controls that does not account for
    class imbalance.

  - "\\integer\\" :  
    `N = \<integer\>`  
    If method is a positive integer, it will be used as N for every row.

- force_new:

  If "Neff" (or "N") already exists in `sumstats_dt`, replace it with
  the recomputed version.

- append_method_name:

  should Neff column have an indicator to explain the method that makes
  it., Default is FALSE unless multiple methods are passed

## Value

A data.table with a new column "Neff" or "N"

## Details

There are many different formulas for calculating ESS, but LDSC is
probably the best method available here, as it doesn't assume that the
proportion of controls:cases is 2:1 (as in GIANT) or 4:1 (as in METAL).
