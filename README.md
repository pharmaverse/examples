# Pharmaverse Examples
This book contains end-to-end examples of using pharmaverse packages together to
achieve common clinical reporting analyses.
The examples use consistent source SDTMs to create ADaMs (such as ADSL and ADAE)
and using these as input to produce some familiar Tables/Listings/Graphs and
associated interactive displays (via Shiny).

Other examples may be included here, e.g. for PK/PD or Therapeutic Area
specifics (such as Oncology or Vaccines).

Note that this examples book should only be used to show how collections of
packages can be used in conjunction - more thorough examples of individual
package usages would always be covered in the package site vignettes and no
need to repeat here.

# Posit Cloud
Each example can be explored via a live and interactive Posit Cloud environment (preconfigured with all required package installations).
Please click the "Launch Posit Cloud" button to work with the live code. <<TO BE ADDED>>
Feel free to copy the repo and customize any examples for your own internal clinical reporting workflows!

# Contributing
Anyone can add examples to this repo for usage of any pharmaverse packages.
To do this make an issue with your idea and the maintainer team will contact
you to provide you the necessary repo access.
One note to dependencies is we don't use `renv` here to manage this examples
book. If you ever need to use a new package the requirement is that you will edit
the root level `DESCRIPTION` file and add it to the Imports. This is a
simple dependency management approach, as then during rendering/deployment our
pipeline will pull the latest released version from CRAN.