## Test environments
* Windows 7, R 3.3.0
* Fedora 20, linux-gnu x86_64, R 3.2.0
* win-builder (devel and release)

## R CMD check results
There was no ERROR or WARNING

There was 1 NOTE:

* Possibly mis-spelled words in DESCRIPTION:  PCA, bayesian, frequentist, logit, morphometric - this are not misspelled.



## Resubmission

> Thanks, we see:

> Possibly mis-spelled words in DESCRIPTION: bayesian (8:461) logit (8:428)
> Can you pls write Bayesian and logistic?

The DESCRIPTION file has been updated.

> More important: All your examples are wrapped in \dontrun{}.Please add examples that actually get executed during the checks so that your code is tested in some way.

It has been corrected.

> Best,
> Uwe Ligges

Thanks for the corrections

Best regards.