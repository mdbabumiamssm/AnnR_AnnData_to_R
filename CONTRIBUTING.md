# Reporting problems and proposing changes

AnnR is authored and maintained by **MD BABU MIA**.

Please open an issue with the AnnR version, `sessionInfo()`,
`annr_build_info()`, operating system, and a minimal synthetic example.
Do not upload patient records or confidential datasets. If possible,
include `h5ad_info()` output after removing identifying names.

For development, install the suggested dependencies, then run:

```r
testthat::test_local(".", stop_on_failure = TRUE)
```

Build and check the source package before proposing a change. See
[validation](docs/VALIDATION.md) for the independent Python comparisons.
Reports are welcome; third-party code contributions must retain accurate
attribution if accepted in a future release.
