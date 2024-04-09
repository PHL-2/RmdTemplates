# RmdTemplates

## Introduction

R markdown templates to initialize as new projects to analyze sequencing runs

## Installation

To install the repo

```
devtools::install_github("PHL-2/RmdTemplates")
```

## Setup

Copy the aux_files directory to a parent directory that will hold all the relevant RProject folders<br/>
Fill in the config_variables.R file in aux_files/config with the appropriate settings

## Creating a template

To use a template for a new run, first generate a new RProject<br/>
Under the new RProject, type the following command on the console

```
rmarkdown::draft(file = "[YYYY-MM-DD]_[run_name]_QC_Report.Rmd", template = "covid_dragen_bs", package = "RmdTemplates", create_dir = FALSE)
```

where [YYYY-MM-DD] is the date and [run_name] is the name of the sequencing run/samples

## Running the analysis

Depending on the template used, add in additional metadata files to the appropriate metadata folders<br/>
Then source the necessary Rscripts to generate additional metadata, submit jobs to process data, and to download the resulting files (this can also be done by running the relevant chunk in [YYYY-MM-DD]\_[run_name]\_QC\_Report.Rmd file)<br/>
After downloading the necessary file, run the 'generate report' chunk to analyze the results
