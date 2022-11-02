# RmdTemplates

## Introduction

R markdown templates to initialize as new projects to analyze sequencing runs

## Installation

To install the repo

```
devtools::install_github("PHL-2/RmdTemplates")
```

## Creating a template

To generate a new Rmd template, type the following command on the console

```
rmarkdown::draft(file = "[YYYY-MM-DD]_[run_name]_QC_Report.Rmd", template = "covid_dragen_bs", package = "RmdTemplates", create_dir = FALSE)
```

where [YYYY-MM-DD] is the date and [run_name] is the name of the sequencing run/samples

## Running the analysis

To create a sequencing report, first fill out the metadata sheet located in the metadata/munge folder<br/>
Then source the generate_barcodes_IDT.R Rscript (can also be done by running the relevant chunk in [YYYY-MM-DD]\_[run_name]\_QC\_Report.Rmd file)<br/>
Then run the last chunk in the [YYYY-MM-DD]\_[run_name]\_QC\_Report.Rmd file
