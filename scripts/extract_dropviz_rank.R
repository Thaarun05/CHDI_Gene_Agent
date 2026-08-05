#!/usr/bin/env Rscript
# extract_dropviz_rank.R
#
# Extract clusters.top from a DropViz rank.zip / rank.Rdata export.
# Writes:
#   clusters_top_raw.csv      — source order and source column names
#   clusters_top_ranked.csv   — descending by target.sum.per.100k
#   top_clusters.json         — compact derived representation
#   rank_extraction_inventory.json
#
# Usage:
#   Rscript extract_dropviz_rank.R <zip_or_rdata> <output_dir>
#
# Exit codes:
#   0 success
#   2 missing clusters.top / unusable input
#   3 validation failure
#   1 other error

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript extract_dropviz_rank.R <zip_or_rdata> <output_dir>\n", file = stderr())
  quit(status = 1)
}

input_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_rdata <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env
}

find_rdata_in_zip <- function(zip_path) {
  members <- unzip(zip_path, list = TRUE)$Name
  rdata <- members[grepl("\\.Rdata$|\\.rda$|\\.RData$", members, ignore.case = TRUE)]
  if (length(rdata) == 0) {
    return(NULL)
  }
  # Prefer rank.Rdata when present.
  prefer <- rdata[grepl("rank\\.Rdata$", rdata, ignore.case = TRUE)]
  if (length(prefer) > 0) prefer[[1]] else rdata[[1]]
}

resolve_input <- function(path) {
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    member <- find_rdata_in_zip(path)
    if (is.null(member)) {
      stop("No .Rdata member found in ZIP: ", path)
    }
    tmp <- tempfile(pattern = "dropviz_rank_", fileext = ".Rdata")
    unzip(path, files = member, junkpaths = TRUE, exdir = dirname(tmp))
    extracted <- file.path(dirname(tmp), basename(member))
    if (!file.exists(extracted)) {
      # junkpaths may have placed it as tmp itself if names collide; try member basename in exdir
      extracted <- file.path(dirname(tmp), basename(member))
    }
    # Ensure stable path
    if (extracted != tmp && file.exists(extracted)) {
      file.rename(extracted, tmp)
    }
    list(rdata = tmp, cleanup = TRUE, member = member)
  } else {
    list(rdata = path, cleanup = FALSE, member = NA_character_)
  }
}

write_json_simple <- function(obj, path) {
  # Minimal JSON writer avoiding package deps.
  esc <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", x)
    x <- gsub("\"", "\\\\\"", x)
    x <- gsub("\n", "\\\\n", x)
    x
  }
  to_json <- function(x) {
    if (is.null(x)) return("null")
    if (is.logical(x) && length(x) == 1) return(if (x) "true" else "false")
    if (is.numeric(x) && length(x) == 1) {
      if (!is.finite(x)) return("null")
      return(as.character(x))
    }
    if (is.character(x) && length(x) == 1) return(paste0("\"", esc(x), "\""))
    if (is.list(x) && is.null(names(x))) {
      return(paste0("[", paste(vapply(x, to_json, character(1)), collapse = ","), "]"))
    }
    if (is.list(x) || is.data.frame(x)) {
      nms <- names(x)
      parts <- vapply(seq_along(x), function(i) {
        paste0("\"", esc(nms[[i]]), "\":", to_json(x[[i]]))
      }, character(1))
      return(paste0("{", paste(parts, collapse = ","), "}"))
    }
    if (is.vector(x)) {
      return(paste0("[", paste(vapply(as.list(x), to_json, character(1)), collapse = ","), "]"))
    }
    paste0("\"", esc(as.character(x)), "\"")
  }
  writeLines(to_json(obj), con = path)
}

ok <- FALSE
inventory <- list(
  input = input_path,
  status = "pending",
  objects = character(0),
  has_clusters_top = FALSE
)

tryCatch({
  resolved <- resolve_input(input_path)
  env <- load_rdata(resolved$rdata)
  objs <- ls(envir = env)
  inventory$objects <- objs
  inventory$rdata_member <- resolved$member

  if (!("clusters.top" %in% objs)) {
    inventory$status <- "missing_clusters_top"
    write_json_simple(inventory, file.path(output_dir, "rank_extraction_inventory.json"))
    cat("clusters.top not found in Rdata. Objects: ", paste(objs, collapse = ", "), "\n", file = stderr())
    quit(status = 2)
  }

  clusters.top <- get("clusters.top", envir = env)
  if (!is.data.frame(clusters.top)) {
    clusters.top <- as.data.frame(clusters.top, stringsAsFactors = FALSE)
  }
  inventory$has_clusters_top <- TRUE
  inventory$nrow <- nrow(clusters.top)
  inventory$ncol <- ncol(clusters.top)
  inventory$columns <- names(clusters.top)

  raw_path <- file.path(output_dir, "clusters_top_raw.csv")
  # Preserve source column names and row order exactly.
  write.csv(clusters.top, file = raw_path, row.names = FALSE)

  if (!("target.sum.per.100k" %in% names(clusters.top))) {
    inventory$status <- "missing_expression_column"
    write_json_simple(inventory, file.path(output_dir, "rank_extraction_inventory.json"))
    cat("target.sum.per.100k missing from clusters.top\n", file = stderr())
    quit(status = 2)
  }

  # Validate CI bounds / finite / nonnegative / nonempty labels
  label_col <- NULL
  for (cand in c("cx.disp", "cluster.disp", "subcluster.disp", "label")) {
    if (cand %in% names(clusters.top)) {
      label_col <- cand
      break
    }
  }
  if (is.null(label_col)) {
    inventory$status <- "rank_validation_failed"
    write_json_simple(inventory, file.path(output_dir, "rank_extraction_inventory.json"))
    cat("No nonempty label column found\n", file = stderr())
    quit(status = 3)
  }

  est <- clusters.top[["target.sum.per.100k"]]
  lo <- if ("target.sum.L.per.100k" %in% names(clusters.top)) clusters.top[["target.sum.L.per.100k"]] else NULL
  hi <- if ("target.sum.R.per.100k" %in% names(clusters.top)) clusters.top[["target.sum.R.per.100k"]] else NULL
  labels <- as.character(clusters.top[[label_col]])

  bad <- which(
    is.na(labels) | !nzchar(trimws(labels)) |
      !is.finite(est) | est < 0 |
      (!is.null(lo) & (!is.finite(lo) | lo < 0)) |
      (!is.null(hi) & (!is.finite(hi) | hi < 0)) |
      (!is.null(lo) & !is.null(hi) & (lo > est | est > hi))
  )
  if (length(bad) > 0 && length(bad) == nrow(clusters.top)) {
    inventory$status <- "rank_validation_failed"
    inventory$invalid_row_count <- length(bad)
    write_json_simple(inventory, file.path(output_dir, "rank_extraction_inventory.json"))
    cat("All rank rows failed validation\n", file = stderr())
    quit(status = 3)
  }

  keep <- setdiff(seq_len(nrow(clusters.top)), bad)
  ranked <- clusters.top[keep, , drop = FALSE]
  ranked <- ranked[order(-ranked[["target.sum.per.100k"]]), , drop = FALSE]
  ranked_path <- file.path(output_dir, "clusters_top_ranked.csv")
  write.csv(ranked, file = ranked_path, row.names = FALSE)

  top_list <- lapply(seq_len(nrow(ranked)), function(i) {
    row <- ranked[i, , drop = FALSE]
    list(
      label = as.character(row[[label_col]][[1]]),
      `target.sum.per.100k` = as.numeric(row[["target.sum.per.100k"]][[1]]),
      `target.sum.L.per.100k` = if (!is.null(lo)) as.numeric(row[["target.sum.L.per.100k"]][[1]]) else NULL,
      `target.sum.R.per.100k` = if (!is.null(hi)) as.numeric(row[["target.sum.R.per.100k"]][[1]]) else NULL,
      gene = if ("gene" %in% names(row)) as.character(row[["gene"]][[1]]) else NULL,
      region = if ("region.disp" %in% names(row)) as.character(row[["region.disp"]][[1]]) else NULL
    )
  })
  top_obj <- list(
    clusters = top_list,
    sort_policy = "descending by target.sum.per.100k; raw source row order preserved in clusters_top_raw.csv; CI validated lower_bound <= estimate <= upper_bound"
  )
  write_json_simple(top_obj, file.path(output_dir, "top_clusters.json"))

  inventory$status <- "success"
  inventory$invalid_row_count <- length(bad)
  inventory$valid_row_count <- length(keep)
  inventory$sort_policy <- top_obj$sort_policy
  write_json_simple(inventory, file.path(output_dir, "rank_extraction_inventory.json"))

  if (resolved$cleanup && file.exists(resolved$rdata)) {
    unlink(resolved$rdata)
  }
  ok <- TRUE
}, error = function(e) {
  cat("extract_dropviz_rank.R error: ", conditionMessage(e), "\n", file = stderr())
  inventory$status <<- "extraction_failed"
  inventory$error <<- conditionMessage(e)
  try(write_json_simple(inventory, file.path(output_dir, "rank_extraction_inventory.json")), silent = TRUE)
  quit(status = 1)
})

if (!ok) quit(status = 1)
quit(status = 0)
