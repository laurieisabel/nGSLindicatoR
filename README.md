nGSLindicatoR: Generate Indicators for Northern Gulf of St. Lawrence
Offshore Communities
================

# Download package

``` r
devtools::install_github("laurieisabel/nGSLindicatoR")
```

# How it works

The package contains two functions:

- `generate_data`: calculate the annual values for each indicator and
  generate the files or data.
- `generate_figures`: generate the figures published in Isabel et
  al. (2025).

``` r
library(nGSLindicatoR)

generate_data(data.directory = "//dcqcimlna01a//BD_Peches//Releves_Poissons_de_Fond_et_Crevette//Donnees_PACES//", repository = "data/", overwrite = TRUE)

generate_figures(repository = "data/", figures = "fig/", overwrite = TRUE)
```

    ## Warning: No shared levels found between `names(values)` of the manual scale and the
    ## data's fill values.
    ## No shared levels found between `names(values)` of the manual scale and the
    ## data's fill values.
