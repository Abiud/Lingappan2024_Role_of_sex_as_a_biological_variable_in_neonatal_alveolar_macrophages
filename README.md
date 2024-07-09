<a target="_blank" href="https://www.lingappanlab.com/" rel="noopener noreferrer"><img align="right" alt="Lingappan Lab logo" src="https://www.lingappanlab.com/_next/static/media/logo.80977683.svg" width="128"/></a>

# Title
Link to article.

The book folder contains quarto files with a step by step analysis.
The R folder conatins all of the functions used on each step.

## GEO

[GEOID](LinkToGEO)

## Docker

The Docker BioContainer Image was used to run this analysis and can be find here: [bioconductor_docker:RELEASE_3_18-R-4.3.2](https://hub.docker.com/layers/bioconductor/bioconductor_docker/RELEASE_3_18-R-4.3.2/images/sha256-6592b272e2cb15ac438dd87154a01d718ea552e27713d4d36315d0b3ae33b8f3?context=explore)

## How to run

R 4.3.2 and the package renv is recommended so the user can accurately replicate the resulsts in the publication.

Clone the project, download the raw counts from GEO and place on raw_data subfolder.

Install needed packages to run pipeline.

```R
  renv::restore()
```

## Questions & Comments

Please contact Abiud (the owner of this repo) for any questions or comments.
