fixture <- function(name) test_path("fixtures", name)

expected <- function() {
  jsonlite::fromJSON(fixture("expected.json"), simplifyVector = TRUE, simplifyMatrix = TRUE)
}

tmp_h5ad <- function() tempfile(fileext = ".h5ad")

# The fixture holds int64 values beyond 2^53 on purpose; reading it warns.
read_fixture <- function(...) {
  withCallingHandlers(read_h5ad(...), warning = function(w) {
    if (grepl("beyond 2\\^53", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}
