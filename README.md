This repository contains the code, text, and results for the manuscript describing the [`hesim`](https://hesim-dev.github.io/hesim/) `R` package. A replication file titled `manuscript.R` can be obtained with: 

```r
knitr::purl(input = "manuscript.rnw", output = "manuscript.R")
```

Package dependencies are manged with the [`renv`](https://rstudio.github.io/renv/articles/renv.html) package. You can install all `R` packages used in the replication file with:

```r
renv::restore()
```

The PDF version of the manuscript can be compiled by running:

```{r}
knitr::knit2pdf(input = "manuscript.rnw", output = "manuscript.tex",
                compiler = "pdflatex")
```

You may alternatively use the "Compile PDF" button within RStudio after modifying *Project Options -> Sweave* so that you weave PDF files with knitr and typeset LaTeX into PDF with pdfLaTeX (you may also modify the *Global Options*).