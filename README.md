# wrapFA: A Wrapper for Factor Analysis using lavaan and MplusAutomation
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/wrapFA?color=brightgreen)](https://cran.r-project.org/package=wrapFA)

## How to cite this package
Nájera, P., Abad, F. J., & Sorrel, M. A. (2022). *wrapFA: A Wrapper for Factor Analysis using lavaan and MplusAutomation*. R package version 0.0.1. https://github.com/Pablo-Najera/wrapFA.
## Features of the package
* Several factor analytic techniques (e.g., CFA, EFA, BSEM) with a simple and unified syntax.
* Determinacy of factor score estimates.
* Bootstrap simulation to evaluate the stability of factor analysis solutions.
## Installation
To install this package from source:
1. Windows users may need to install the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and include the checkbox option of installing Rtools to their path for easier command line usage. Mac users will have to download the necessary tools from the [Xcode](https://apps.apple.com/ca/app/xcode/id497799835?mt=12) and its related command line tools (found within Xcode's Preference Pane under Downloads/Components); most Linux distributions should already have up to date compilers (or if not they can be updated easily).
2. Install the `devtools` package (if necessary), and install the package from the Github source code.

```
#install.packages("devtools")
devtools::install_github("Pablo-Najera/wrapFA")
```

## Bug reports
Please report any bugs at https://github.com/Pablo-Najera/wrapFA/issues.
