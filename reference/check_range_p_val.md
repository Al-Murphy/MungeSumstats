# Ensure that the p values are not \>1 and if so set to 1

Ensure that the p values are not \>1 and if so set to 1

## Usage

``` r
check_range_p_val(sumstats_dt, convert_large_p, convert_neg_p, imputation_ind)
```

## Source

` sumstats_dt <- MungeSumstats:::formatted_example() sumstats_dt$P[1:3] <- 5 sumstats_dt$P[6:10] <- -5 sumstats <- check_range_p_val(sumstats_dt = sumstats_dt, convert_large_p = TRUE, convert_neg_p = TRUE, imputation_ind = TRUE) `

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS

- convert_large_p:

  Binary, should p-values \>1 be converted to 1? P-values \>1 should not
  be possible and can cause errors with LDSC/MAGMA and should be
  converted. Default is TRUE.

- convert_neg_p:

  Binary, should p-values \<0 be converted to 0? Negative p-values
  should not be possible and can cause errors with LDSC/MAGMA and should
  be converted. Default is TRUE.

- imputation_ind:

  Binary Should a column be added for each imputation step to show what
  SNPs have imputed values for differing fields. This includes a field
  denoting SNP allele flipping (flipped). On the flipped value, this
  denoted whether the alelles where switched based on MungeSumstats
  initial choice of A1, A2 from the input column headers and thus may
  not align with what the creator intended.**Note** these columns will
  be in the formatted summary statistics returned. Default is FALSE.

## Value

list containing sumstats_dt, the modified summary statistics data table
object
