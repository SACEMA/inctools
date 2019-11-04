# Incidence Estimation Tools (*inctools*)

Tools for estimating HIV incidence from biomarker data, and calibrating recent infection tests.

This package is a more robust implementation of the deprecated spreadsheet tools released as *ABIE* and also contains the functionality previously available in the deprecated package *ritcalib*.

*inctools* is available on the [Comprehensive R Archive Network](https://cran.r-project.org/web/packages/inctools/index.html). You can also find the latest stable release under [Releases](https://github.com/SACEMA/inctools/releases). Releases on this repository may be slightly ahead of the version available through CRAN, but we recommend installing from CRAN.

Older spreadsheet-based tools for inferring incidence using biomarkers for 'recent infection' can be found [here](http://www.incidence-estimation.org/page/tools-for-incidence-from-biomarkers-for-recent-infection). These tools are largely superceded by *inctools*.

## About

These tools are released as an R package by the [DST/NRF Centre of Excellence in Epidemiological Modelling and Analysis (SACEMA)](http://www.sacema.org/) in collaboration with the [Consortium for the Evaluation and Performance of HIV Incidence Assays (CEPHIA)](http://www.incidence-estimation.org/page/cephia) and other partners. The package maintainer can be contacted at [eduardgrebe@sun.ac.za](mailto:eduardgrebe@sun.ac.za).

## Installation

Install released version from CRAN as follows:

```
install.packages("inctools")
```

### Current development version

To build and install the current development version, make sure you have
`devtools` installed and execute the following line:

```
devtools::install_github("eduardgrebe/inctools/inctools")
```

## Usage terms

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Contributors

Several people contributed to this package. We encourage contributions from the
community of users.

Developers include Alex Welte, Eduard Grebe, Avery McIntosh, Petra Bäumler,
Simon Daniel and Yuruo Li, with contributions by Cari van Schalkwyk,
Reshma Kassanjee, Hilmarie Brand, Stefano Ongarello and Yusuke Asai.

## Contact

For usage support and bug reports, please contact the package maintainer,
Eduard Grebe, or submit issues on this repository.

## References

For more information on incidence inference based on cross-sectionally obtained biomarkers for 'recent infection', see [this page](http://www.incidence-estimation.org/page/theory-cross-sectional-recent-infection-tests), and [this timeline](http://www.incidence-estimation.org/page/timeline) with links to references, as well as [a brief review with references](http://www.incidence-estimation.org/page/theory-review-and-references-incidence-inference-using-biomarkers-for-recent-infection). The key reference for biomarker-based incidence estimation is [Kassanjee, McWalter, Bärnighausen & Welte (2012)](http://dx.doi.org/10.1097/EDE.0b013e3182576c07).
