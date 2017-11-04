# DGSA

This package implements distance based generalized sensitivity analysis method (DGSA).
For more information about the method please read:

>Fenwick, Darryl, Céline Scheidt, and Jef Caers. "Quantifying asymmetric parameter
interactions in sensitivity analysis: Application to reservoir modeling." Mathematical
Geosciences 46.4 (2014): 493-511.

You can install the package manually or from GitHub:

## Manual installation

Clone to your own computer, open Rstudio or R and type:

```R
install.packages("package_Location", repos = NULL, type=source)
```

You also have an option to load your clone as a project in Rstudio and then proceed with
project build <kbd>Ctrl</kbd>+<kbd>Shift</kbd>+<kbd>B</kbd>. This is particularly useful
if you plan on adding/modifying features.

## Installing from GitHub

First you will need to make sure that you have `devtools` installed, then you can install
with the following command:

```R
devtools::install_github("ogru/DGSA")
```

After that you can proceed to load the routines in the usual way with `library(DGSA)`.

## Examples of usage

To view some of the package capabilities go to the following link:

https://rawgit.com/ogru/DGSA/master/vignette/Intro_to_DGSA.html

It will view the vignette without downloading the package/file.

## Related software

If you prefer to work with MATLAB, you can find the original MATLAB implementation here:
https://github.com/scheidtc/dGSA. However, without the matrix plot for visualization of
all sensitivities/importances.

## Citation

If you find this sofware useful please consider citing it in your publications. You may use the following bibtex type citation entry:
```
@Manual{DGSA_2016,
  title  = {DGSA, an \proglang{R} package},
  author = {Ognjen Grujic},
  year   = {2016},
  note   = {\proglang{R}~package version~1.0},
  url    = {https://github.com/ogru/DGSA},
}
```
