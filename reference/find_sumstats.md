# Search Open GWAS for datasets matching criteria

For each argument, searches for any datasets matching a case-insensitive
substring search in the respective metadata column. Users can supply a
single character string or a list/vector of character strings.

## Usage

``` r
find_sumstats(
  ids = NULL,
  traits = NULL,
  years = NULL,
  consortia = NULL,
  authors = NULL,
  populations = NULL,
  categories = NULL,
  subcategories = NULL,
  builds = NULL,
  pmids = NULL,
  min_sample_size = NULL,
  min_ncase = NULL,
  min_ncontrol = NULL,
  min_nsnp = NULL,
  include_NAs = FALSE
)
```

## Arguments

- ids:

  List of Open GWAS study IDs (e.g. `c("prot-a-664", "ieu-b-4760")`).

- traits:

  List of traits (e.g. `c("parkinson", "Alzheimer")`).

- years:

  List of years (e.g. `seq(2015,2021)` or `c(2010, 2012, 2021)`).

- consortia:

  List of consortia (e.g. `c("MRC-IEU","Neale Lab")`.

- authors:

  List of authors (e.g. `c("Elsworth","Kunkle","Neale")`).

- populations:

  List of populations (e.g. `c("European","Asian")`).

- categories:

  List of categories (e.g.
  `c("Binary","Continuous","Disease","Risk factor"))`).

- subcategories:

  List of categories (e.g. `c("neurological","Immune","cardio"))`).

- builds:

  List of genome builds (e.g. `c("hg19","grch37")`).

- pmids:

  List of PubMed ID (exact matches only) (e.g.
  `c(29875488, 30305740, 28240269)`).

- min_sample_size:

  Minimum total number of study participants (e.g. `5000`).

- min_ncase:

  Minimum number of case participants (e.g. `1000`).

- min_ncontrol:

  Minimum number of control participants (e.g. `1000`).

- min_nsnp:

  Minimum number of SNPs (e.g. `200000`).

- include_NAs:

  Include datasets with missing metadata for size criteria (i.e.
  `min_sample_size`, `min_ncase`, or `min_ncontrol`).

## Value

(Filtered) GWAS metadata table.

## Details

To authenticate, you need to generate a token from the OpenGWAS website.
The token behaves like a password, and it will be used to authorise the
requests you make to the OpenGWAS API. Here are the steps to generate
the token and then have `ieugwasr` automatically use it for your
queries:

1.  Login to https://api.opengwas.io/profile/

2.  Generate a new token

3.  Add `OPENGWAS_JWT=<token>` to your .Renviron file, thi can be edited
    in R by running
    [`usethis::edit_r_environ()`](https://usethis.r-lib.org/reference/edit.html)

4.  Restart your R session

5.  To check that your token is being recognised, run
    [`ieugwasr::get_opengwas_jwt()`](https://mrcieu.github.io/ieugwasr/reference/get_opengwas_jwt.html).
    If it returns a long random string then you are authenticated.

6.  To check that your token is working, run
    [`ieugwasr::user()`](https://mrcieu.github.io/ieugwasr/reference/user.html).
    It will make a a request to the API for your user information using
    your token. It should return a list with your user information. If
    it returns an error, then your token is not working.

7.  Make sure you have submitted use

By default, returns metadata for all studies currently in Open GWAS
database.

## Examples

``` r
# Only run the examples if user has internet access
# and if access token has been added
if(try(is.character(getURL("www.google.com")))==TRUE && ieugwasr::get_opengwas_jwt()!=""){
### By ID
metagwas <- find_sumstats(ids = c(
    "ieu-b-4760",
    "prot-a-1725",
    "prot-a-664"
))
### By ID and sample size
metagwas <- find_sumstats(
    ids = c("ieu-b-4760", "prot-a-1725", "prot-a-664"),
    min_sample_size = 5000
)
### By criteria
metagwas <- find_sumstats(
    traits = c("alzheimer", "parkinson"),
    years = seq(2015, 2021)
)
}
#> Error in getURL("www.google.com") : could not find function "getURL"
```
