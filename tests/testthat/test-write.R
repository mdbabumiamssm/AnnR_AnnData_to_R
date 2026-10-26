make_annr <- function() {
  set.seed(1)
  n_obs <- 9
  n_var <- 6
  X <- Matrix::rsparsematrix(n_obs, n_var, density = 0.4)
  X <- methods::as(X, "RsparseMatrix")
  obs <- data.frame(
    group = factor(sample(c("a", "b"), n_obs, TRUE)),
    level = factor(sample(c("lo", "hi"), n_obs, TRUE), levels = c("lo", "hi"), ordered = TRUE),
    score = rnorm(n_obs),
    count = sample.int(100, n_obs),
    flag = rep(c(TRUE, FALSE, NA), 3),
    name = c(letters[1:8], NA),
    maybe = c(1L, NA, 3L, 4L, 5L, NA, 7L, 8L, 9L),
    row.names = paste0("cell", seq_len(n_obs))
  )
  var <- data.frame(hv = rep(c(TRUE, FALSE), 3), row.names = paste0("gene", seq_len(n_var)))
  AnnR(
    X = X, obs = obs, var = var,
    layers = list(dense = matrix(rnorm(n_obs * n_var), n_obs), ints = matrix(1:(n_obs * n_var), n_obs)),
    obsm = list(X_pca = matrix(rnorm(n_obs * 3), n_obs)),
    varm = list(PCs = matrix(rnorm(n_var * 3), n_var)),
    obsp = list(knn = Matrix::rsparsematrix(n_obs, n_obs, 0.2)),
    uns = list(title = "t", k = 15L, x = 2.5, flag = TRUE, colors = I("#ffffff"),
               nested = list(a = 1:3, b = list(c = "deep")), cplx = complex(real = 1, imaginary = -2),
               arr3 = array(1:24, c(2, 3, 4)), strmat = matrix(c("a", "b", "c", "d"), 2))
  )
}

expect_same_annr <- function(a, b) {
  attr(a, "h5ad_types") <- NULL
  attr(b, "h5ad_types") <- NULL
  expect_equal(as.matrix(a$X), as.matrix(b$X), ignore_attr = TRUE)
  expect_identical(a$obs, b$obs)
  expect_identical(a$var, b$var)
  expect_equal(a$layers, b$layers, ignore_attr = TRUE)
  expect_equal(a$obsm, b$obsm, ignore_attr = TRUE)
  expect_equal(a$varm, b$varm, ignore_attr = TRUE)
  expect_equal(lapply(a$obsp, as.matrix), lapply(b$obsp, as.matrix), ignore_attr = TRUE)
}

for (comp in c("none", "gzip", "lzf", "zstd")) {
  test_that(paste("R objects round-trip with compression", comp), {
    if (comp == "zstd") skip_if(annr_build_info()[["zstd"]] != "yes", "built without zstd")
    ad <- make_annr()
    f <- tmp_h5ad()
    write_h5ad(ad, f, compression = comp)
    back <- read_h5ad(f)
    expect_same_annr(ad, back)
    expect_identical(back$uns$title, "t")
    expect_identical(back$uns$k, 15L)
    expect_identical(back$uns$x, 2.5)
    expect_identical(back$uns$flag, TRUE)
    expect_identical(back$uns$colors, "#ffffff")
    expect_identical(back$uns$nested$a, 1:3)
    expect_identical(back$uns$nested$b$c, "deep")
    expect_identical(back$uns$cplx, complex(real = 1, imaginary = -2))
    expect_identical(back$uns$arr3, array(1:24, c(2, 3, 4)))
    expect_identical(back$uns$strmat, matrix(c("a", "b", "c", "d"), 2))
    info <- h5ad_info(f)
    expect_identical(info$encoding[info$path == "obs/maybe"], "nullable-integer")
    expect_identical(info$encoding[info$path == "obs/flag"], "nullable-boolean")
    expect_identical(info$encoding[info$path == "obs/name"], "nullable-string-array")
    expect_identical(info$encoding[info$path == "uns/k"], "numeric-scalar")
    expect_identical(info$encoding[info$path == "uns/colors"], "string-array")
    expect_identical(info$dtype[info$path == "obs/group/codes"], "int8")
  })
}

test_that("re-writing a file reproduces its dtypes", {
  ad <- read_fixture(fixture("modern.h5ad"))
  f <- tmp_h5ad()
  write_h5ad(ad, f)
  a <- h5ad_info(fixture("modern.h5ad"))
  b <- h5ad_info(f)
  common <- intersect(a$path[a$kind == "dataset"], b$path[b$kind == "dataset"])
  expect_gt(length(common), 30)
  expect_identical(b$dtype[match(common, b$path)], a$dtype[match(common, a$path)])
  expect_identical(b$encoding[match(common, b$path)], a$encoding[match(common, a$path)])
})

test_that("a failed write leaves no partial file behind", {
  ad <- make_annr()
  ad$uns$bad <- list(`a/b` = 1)
  f <- tmp_h5ad()
  expect_error(write_h5ad(ad, f), "may not contain")
  expect_false(file.exists(f))
  expect_length(list.files(dirname(f), pattern = "^\\.annr_"), 0)
})

test_that("subsetting keeps every slot aligned", {
  ad <- make_annr()
  s <- ad[c(2, 4, 6), c("gene1", "gene3")]
  expect_equal(dim(s), c(3L, 2L))
  expect_identical(obs_names(s), paste0("cell", c(2, 4, 6)))
  expect_equal(dim(s$layers$dense), c(3L, 2L))
  expect_equal(nrow(s$obsm$X_pca), 3L)
  expect_equal(dim(s$obsp$knn), c(3L, 3L))
  expect_equal(nrow(s$varm$PCs), 2L)
  f <- tmp_h5ad()
  write_h5ad(s, f)
  expect_equal(dim(read_h5ad(f)), c(3L, 2L))
})
