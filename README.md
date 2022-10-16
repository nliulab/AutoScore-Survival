AutoScore-Survival: Developing interpretable machine learning-based
time-to-event scores with right-censored survival data
================

AutoScore-Survival is an method extension to AutoScore,  and a novel machine learning framework to automate the
development of interpretable time-to-event scores. AutoScore-Survival
consists of six modules: 1) variable ranking with machine learning, 2)
variable transformation, 3) score derivation, 4) model selection, 5)
domain knowledge-based score fine-tuning, and 6) performance evaluation.
AutoScore-Survival could seamlessly
generate risk scores based on survival data, which can be easily
implemented and validated in clinical practice. Moreover, it enables
users to build transparent and interpretable time-to-event scores
quickly in a straightforward manner.

> ### AutoScore-Survival has been merged with the [AutoScore package](https://github.com/nliulab/AutoScore/). Please visit [AutoScore bookdown page](https://nliulab.github.io/AutoScore/) for a full tutorial.

## Citation

Xie F, Ning Y, Yuan H, Goldstein BA, Ong MEH, Liu N, Chakraborty B. 
AutoScore-Survival: Developing interpretable machine learning-based time-to-event scores with right-censored survival data
Journal of biomedical informatics, 125 (2022) (<https://doi.org/10.1016/j.jbi.2021.103959>)

Xie F, Chakraborty B, Ong MEH, Goldstein BA, Liu N. AutoScore: A Machine
Learning-Based Automatic Clinical Score Generator and Its Application to
Mortality Prediction Using Electronic Health Records. JMIR Medical
Informatics 2020;8(10):e21798 (<http://dx.doi.org/10.2196/21798>)

## Contact

  - Feng Xie (Email: <xief@u.duke.nus.edu>)
  - Nan Liu (Email: <liu.nan@duke-nus.edu.sg>)

## Package installation

Install from GitHub or CRAN:

``` r
# From Github
install.packages("devtools")
library(devtools)
install_github(repo = "nliulab/AutoScore", build_vignettes = TRUE)

# From CRAN (recommended)
install.packages("AutoScore")
```

Load AutoScore package: 

``` r
library(AutoScore)
```



