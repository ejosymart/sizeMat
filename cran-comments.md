## Test environments
* Windows 7, R 3.3.0
* Fedora 20, linux-gnu x86_64, R 3.2.0
* win-builder (devel and release)

## R CMD check results
There was no ERROR or WARNING

There was 1 NOTE:

* Possibly mis-spelled words in DESCRIPTION:  PCA, bayesian, frequentist, logit, morphometric - this are not misspelled.



## Resubmission (first round)

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



## Resubmission (second round)

> Now we're seeing

> ** running examples for arch 'i386' ... [43s] NOTE
> Examples with CPU or elapsed time > 10s
>                 user system elapsed
> print.gonadMat 10.70   0.03   10.74
> gonad_mature   10.50   0.05   10.54
> plot.gonadMat  10.52   0.02   10.56
> ** running examples for arch 'x64' ... [42s] NOTE
> Examples with CPU or elapsed time > 10s
>                 user system elapsed
> gonad_mature   10.62   0.08   10.70
> print.gonadMat 10.37   0.02   10.39

> Could you please shorten those examples, so they run in about 5s each?
> Duncan Murdoch


It has been corrected.

Thanks for the corrections.
Best regards.