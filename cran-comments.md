# CRAN comments for fastDist 1.1.0

## Submission

This is the first submission of fastDist to CRAN.

## Test environments

* local Linux (Ubuntu), R release
* win-builder (devel and release)
* R-hub: Windows Server (R devel), Ubuntu Linux (R release), Fedora Linux
  (R devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

(There may be a NOTE about a new submission on win-builder, which is
expected for a first release.)

## Notes for the reviewers

* The package contains C++ code parallelised with OpenMP. OpenMP is enabled
  through the portable `$(SHLIB_OPENMP_CXXFLAGS)` mechanism, so the package
  also builds and runs correctly on toolchains without OpenMP support.
* Examples, tests and the vignette only use small matrices, so check times
  are short.
