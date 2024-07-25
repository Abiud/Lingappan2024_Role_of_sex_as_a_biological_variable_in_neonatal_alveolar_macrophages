<a target="_blank" href="https://www.lingappanlab.com/" rel="noopener noreferrer"><img align="right" alt="Lingappan Lab logo" src="https://www.lingappanlab.com/_next/static/media/logo.80977683.svg" width="128"/></a>

# Impact of neonatal hyperoxia exposure on lung myeloid cells

The book folder contains quarto files with a step by step analysis.
The R folder conatins all of the functions used on each step.

## GEO

[GSE267573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267573)

## Docker

The Docker BioContainer Image was used to run this analysis and can be found here: [bioconductor_docker:3.18-R-4.3.2](https://hub.docker.com/layers/bioconductor/bioconductor_docker/3.18-R-4.3.2/images/sha256-fd5a50d01bdf723396338dc6e9437ad059d50da0da2740e51d7438c8c2787721?context=explore)

## How to run

R 4.3.2 and the package renv is recommended so the user can accurately replicate the resulsts in the publication.

Clone the project, download the raw counts (GSE267573_raw_counts.csv.gz) from GEO gunzip it and place on raw_data subfolder.

> I like to create a new container using the image mentioned above and mounting the folder in my host machine to /home in the container. Then I attach visual studio code (with the R and quarto extensions installed) to the mounted folder (using the visual studio docker extension) and run the pipeline with these steps:

Install needed packages to run pipeline.

In the R console:

```R
install.packages("renv")

renv::init() # Restore the project from the provided lockfile.
```

In the terminal:

```bash
quarto render book/
```

## Questions & Comments

Please contact Abiud (the owner of this repo) for any questions or comments.
