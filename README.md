# mammalitv
Data processing and analysis from the Ecography paper on NEON small mammals

Code associated with Read, Q. D., J. M. Grady, P. L. Zarnetske, S. Record, B. Baiser, J. Belmaker, M.-N. Tuanmu, A. Strecker, L. Beaudrot, and K. M. Thibault. 2018. Among-species overlap in rodent body size distributions predicts species richness along a temperature gradient. *Ecography*. DOI: 10.1111/ecog.03641

## Location of data

The data is mostly stored on MSU's HPCC server. The code can be run remotely on that server.

## Description of subdirectories

### code

- **analysis**: Functions to calculate overlap metrics, which are implemented for NEON's small mammal data.
- **ms_revisions**: Scripts to do analyses specifically requested by manuscript reviewers.
- **vis**: Code to make plots and other data visualizations

### results

This directory has .Rmd files which describe R code and make visualizations. The .Rmd files refer to files stored locally on the user's hard drive so the files will not run as is.


*This readme last modified by QDR, 26 April 2018.*