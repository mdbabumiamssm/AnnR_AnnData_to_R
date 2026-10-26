#' Convert an AnnR object to a Seurat object
#'
#' Builds a Seurat v5 object (`SeuratObject` >= 5.0). The `counts` and
#' `data` layers of the assay are taken from the AnnR matrices named by
#' `counts` and `data` (`"X"` or a layer name; `NULL` to omit); every other
#' layer is kept as an assay layer of the same name. `obs` becomes
#' the cell metadata, `var` the feature metadata, 2-D numeric `obsm` entries
#' become dimensional reductions (`X_pca` -> `pca`, ...), `obsp` entries become
#' graphs and `raw$X` becomes an assay named `"raw"`. Seurat forbids `_` in
#' feature names and duplicated names; such names are adjusted and the
#' originals kept in `misc$AnnR` so that [write_h5ad()] restores them, together
#' with `uns`, `varm` and the h5ad dtype records.
#'
#' @param x An [AnnR] object.
#' @param assay Assay name.
#' @param counts,data Matrices for the `counts` and `data` layers.
#' @param project Project name.
#' @param ... Unused.
#' @return A `Seurat` object.
#' @export
as_Seurat <- function(x, assay = "RNA", counts = "counts", data = "X", project = "AnnR", ...) {
  .require("SeuratObject")
  fm <- isTRUE(attr(x, "annr_feature_major"))
  types <- attr(x, "h5ad_types")
  x <- unclass(x)
  cells <- rownames(x$obs)
  genes <- rownames(x$var)
  pick <- function(name) {
    if (is.null(name)) return(NULL)
    if (identical(name, "X")) return(x$X)
    x$layers[[name]]
  }
  seu_cells <- .seurat_names(cells)
  seu_genes <- .seurat_names(genes, features = TRUE)
  prep <- function(m) {
    m <- .to_feature_major(m, fm)
    if (methods::is(m, "Matrix")) {
      m@Dimnames <- list(seu_genes, seu_cells)
    } else {
      dimnames(m) <- list(seu_genes, seu_cells)
    }
    m
  }
  # Seurat layer name -> AnnR matrix ("X" or a layer key)
  layer_map <- character()
  if (!is.null(pick(data))) layer_map[["data"]] <- data
  if (!is.null(pick(counts)) && !identical(counts, data)) layer_map[["counts"]] <- counts
  if (!is.null(x$X) && !("X" %in% layer_map)) layer_map[["X"]] <- "X"
  for (k in names(x$layers)) {
    if (k %in% layer_map) next
    nm <- if (k %in% names(layer_map)) paste0(k, "_layer") else k
    layer_map[[nm]] <- k
  }
  x_absent <- is.null(x$X)
  # Seurat drops attributes on layers: keep exact 64-bit integer patterns aside
  int64_bits <- list()
  for (an in unique(layer_map)) {
    b <- attr(pick(an), "h5ad_int64")
    if (!is.null(b)) {
      m0 <- pick(an)
      int64_bits[[an]] <- if (fm) b else as.vector(t(matrix(b, nrow(m0), ncol(m0))))
    }
  }
  mats <- lapply(layer_map, function(an) prep(pick(an)))
  if (!length(mats)) {
    # Seurat needs a matrix; record that the AnnR object had none.
    mats <- list(counts = Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                               dims = c(length(genes), length(cells)),
                                               dimnames = list(seu_genes, seu_cells)))
  }
  old <- options(Seurat.object.assay.calcn = FALSE)
  on.exit(options(old), add = TRUE)
  first <- names(mats)[1L]
  a5 <- if (identical(first, "counts")) SeuratObject::CreateAssay5Object(counts = mats[[first]])
        else SeuratObject::CreateAssay5Object(data = mats[[first]])
  if (!identical(first, "counts") && !identical(first, "data")) {
    # CreateAssay5Object stored it as "data"; rename by re-adding
    SeuratObject::LayerData(a5, layer = first) <- mats[[first]]
    SeuratObject::LayerData(a5, layer = "data") <- NULL
  }
  for (nm in names(mats)[-1L]) SeuratObject::LayerData(a5, layer = nm) <- mats[[nm]]
  meta <- x$obs
  rownames(meta) <- seu_cells
  seu <- SeuratObject::CreateSeuratObject(a5, assay = assay, meta.data = meta, project = project)
  if (ncol(x$var)) {
    fmeta <- x$var
    rownames(fmeta) <- seu_genes
    seu[[assay]] <- SeuratObject::AddMetaData(seu[[assay]], metadata = fmeta)
  }
  used <- c(assay, "raw")
  unique_name <- function(nm, suffix) {
    while (nm %in% used) nm <- paste0(nm, suffix)
    used <<- c(used, nm)
    nm
  }
  extra_obsm <- list()
  red_map <- character()
  for (k in names(x$obsm)) {
    m <- x$obsm[[k]]
    if (!is.matrix(m) || !is.numeric(m) || ncol(m) == 0L) {
      extra_obsm[[k]] <- m
      next
    }
    name <- unique_name(sub("^X_", "", k), "_obsm")
    key <- gsub("[^A-Za-z0-9]", "", name)
    if (!nzchar(key) || !grepl("^[A-Za-z]", key)) key <- paste0("X", key)
    key <- paste0(key, "_")
    emb <- m
    dimnames(emb) <- list(seu_cells, paste0(key, seq_len(ncol(emb))))
    loadings <- new(Class = "matrix")
    if (identical(k, "X_pca") && is.matrix(x$varm$PCs) && ncol(x$varm$PCs) == ncol(emb)) {
      loadings <- x$varm$PCs
      dimnames(loadings) <- list(seu_genes, colnames(emb))
    }
    seu[[name]] <- SeuratObject::CreateDimReducObject(embeddings = emb, loadings = loadings,
                                                      key = key, assay = assay)
    red_map[[name]] <- k
  }
  graph_map <- character()
  for (k in names(x$obsp)) {
    name <- unique_name(paste0(assay, "_", k), "_obsp")
    g <- .as_dgC(x$obsp[[k]])
    g@Dimnames <- list(seu_cells, seu_cells)
    seu[[name]] <- SeuratObject::as.Graph(g)
    graph_map[[name]] <- k
  }
  if (!is.null(x$raw) && !is.null(x$raw$X)) {
    rg <- rownames(x$raw$var %||% data.frame(row.names = .default_names(ncol(x$raw$X))))
    rm <- .to_feature_major(x$raw$X, fm)
    if (methods::is(rm, "Matrix")) rm@Dimnames <- list(.seurat_names(rg, TRUE), seu_cells)
    else dimnames(rm) <- list(.seurat_names(rg, TRUE), seu_cells)
    seu[["raw"]] <- SeuratObject::CreateAssay5Object(data = rm)
  }
  seu@misc <- x$uns
  seu@misc$AnnR <- list(
    types = types, obs_names = cells, var_names = genes, reductions = red_map,
    graphs = graph_map, layer_map = layer_map, x_absent = x_absent, int64_bits = int64_bits,
    obsm_other = extra_obsm, varm = x$varm, varp = x$varp, raw = x$raw,
    obs_columns = names(x$obs), assay = assay,
    seurat_cells = colnames(seu), seurat_genes = rownames(seu[[assay]])
  )
  seu
}

.seurat_names <- function(n, features = FALSE) {
  n <- as.character(n)
  n[is.na(n) | !nzchar(n)] <- "unnamed"
  if (features) n <- gsub("_", "-", n, fixed = TRUE)
  make.unique(n)
}

#' @rdname as_AnnR
#' @param assay For `Seurat`: assay to convert (default: the default assay).
#' @param X For `Seurat`: layer used as `X` (default `"data"`).
#' @param layers For `Seurat`: further layers to keep (default `"counts"`).
#' @export
as_AnnR.Seurat <- function(x, assay = NULL, X = "data", layers = "counts", ...) {
  .require("SeuratObject")
  info <- x@misc$AnnR
  assay <- assay %||% info$assay %||% SeuratObject::DefaultAssay(x)
  a <- x[[assay]]
  if (!is.null(info$seurat_cells)) {
    if (!identical(assay, info$assay) ||
        !identical(colnames(x), info$seurat_cells) ||
        !identical(rownames(a), info$seurat_genes)) {
      stop("the converted Seurat object's assay or cell/feature identities have changed; ",
           "subset or rename the AnnR object before converting to Seurat, ",
           "so saved metadata cannot be assigned to the wrong cells or features", call. = FALSE)
    }
  }
  have <- SeuratObject::Layers(a)
  get <- function(l) SeuratObject::LayerData(x, assay = assay, layer = l)
  cells <- colnames(x)
  genes <- rownames(a)
  if (!is.null(info$obs_names) && length(info$obs_names) == length(cells)) cells <- info$obs_names
  if (!is.null(info$var_names) && length(info$var_names) == length(genes)) genes <- info$var_names
  Xm <- NULL
  lay <- list()
  take <- function(l) .strip_dimnames(.to_obs_major(get(l)))
  if (!is.null(info$layer_map)) {
    # invert the mapping recorded by as_Seurat(): Seurat layer -> AnnR matrix
    for (sl in names(info$layer_map)) {
      an <- info$layer_map[[sl]]
      if (!(sl %in% have)) next
      m <- take(sl)
      b <- info$int64_bits[[an]]
      if (!is.null(b) && !methods::is(m, "Matrix")) {
        # bits were stored features x cells; m is cells x features
        attr(m, "h5ad_int64") <- as.vector(t(matrix(b, ncol(m), nrow(m))))
      }
      if (identical(an, "X")) Xm <- m else lay[[an]] <- m
    }
    if (isTRUE(info$x_absent)) Xm <- NULL
  } else {
    if (!is.null(X) && X %in% have) Xm <- take(X)
    for (l in setdiff(intersect(layers, have), X)) {
      if (nrow(get(l)) == length(genes)) lay[[l]] <- take(l)
    }
  }
  obs <- x[[]]
  if (!is.null(info$obs_columns)) {
    keep <- intersect(info$obs_columns, names(obs))
    extra <- setdiff(names(obs), c(info$obs_columns, "orig.ident"))
    obs <- obs[, c(keep, extra), drop = FALSE]
  }
  rownames(obs) <- cells
  fmeta <- tryCatch(a[[]], error = function(e) NULL)
  var <- if (is.data.frame(fmeta) && nrow(fmeta) == length(genes)) fmeta else data.frame(row.names = seq_along(genes))
  rownames(var) <- genes
  obsm <- list()
  for (r in SeuratObject::Reductions(x)) {
    k <- if (!is.null(info$reductions) && r %in% names(info$reductions)) info$reductions[[r]] else paste0("X_", r)
    obsm[[k]] <- .strip_dimnames(SeuratObject::Embeddings(x, reduction = r))
  }
  for (k in names(info$obsm_other)) obsm[[k]] <- info$obsm_other[[k]]
  varm <- info$varm %||% list()
  obsp <- list()
  for (g in SeuratObject::Graphs(x)) {
    k <- if (!is.null(info$graphs) && g %in% names(info$graphs)) info$graphs[[g]] else g
    obsp[[k]] <- .strip_dimnames(methods::as(x[[g]], "dgCMatrix"))
  }
  uns <- x@misc
  uns$AnnR <- NULL
  raw <- info$raw
  out <- AnnR(X = Xm, obs = obs, var = var, obsm = obsm, varm = varm, obsp = obsp,
              varp = info$varp %||% list(), layers = lay, uns = uns, raw = raw)
  attr(out, "h5ad_types") <- info$types
  out
}
