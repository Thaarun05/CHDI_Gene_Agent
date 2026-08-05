#!/usr/bin/env Rscript
# inspect_dropviz_rdata.R
#
# Generic inventory/export helper for DropViz ZIP / Rdata packages (t-SNE etc.).
# Writes:
#   rdata_inventory.json
#   optional CSV dumps for data.frame objects (safe names only)
#
# Usage:
#   Rscript inspect_dropviz_rdata.R <zip_or_rdata> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript inspect_dropviz_rdata.R <zip_or_rdata> <output_dir>\n", file = stderr())
  quit(status = 1)
}

input_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_json_simple <- function(obj, path) {
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
    if (is.list(x)) {
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

find_rdata_in_zip <- function(zip_path) {
  members <- unzip(zip_path, list = TRUE)$Name
  rdata <- members[grepl("\\.Rdata$|\\.rda$|\\.RData$", members, ignore.case = TRUE)]
  if (length(rdata) == 0) return(NULL)
  prefer <- rdata[grepl("tsne\\.Rdata$|rank\\.Rdata$|scatter\\.Rdata$", rdata, ignore.case = TRUE)]
  if (length(prefer) > 0) prefer[[1]] else rdata[[1]]
}

resolve_input <- function(path) {
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    member <- find_rdata_in_zip(path)
    if (is.null(member)) stop("No .Rdata member found in ZIP: ", path)
    tmp <- tempfile(pattern = "dropviz_inspect_", fileext = ".Rdata")
    unzip(path, files = member, junkpaths = TRUE, exdir = dirname(tmp))
    extracted <- file.path(dirname(tmp), basename(member))
    if (file.exists(extracted) && extracted != tmp) file.rename(extracted, tmp)
    list(rdata = tmp, cleanup = TRUE, member = member, zip_members = unzip(path, list = TRUE)$Name)
  } else {
    list(rdata = path, cleanup = FALSE, member = NA_character_, zip_members = character(0))
  }
}

safe_name <- function(name) {
  gsub("[^A-Za-z0-9._-]+", "_", name)
}

tryCatch({
  resolved <- resolve_input(input_path)
  env <- new.env(parent = emptyenv())
  load(resolved$rdata, envir = env)
  objs <- ls(envir = env)

  object_info <- lapply(objs, function(nm) {
    val <- get(nm, envir = env)
    info <- list(
      name = nm,
      class = paste(class(val), collapse = ","),
      typeof = typeof(val)
    )
    if (is.data.frame(val)) {
      info$nrow <- nrow(val)
      info$ncol <- ncol(val)
      info$columns <- names(val)
      # Export a CSV dump for inspection (cap rows).
      out_csv <- file.path(output_dir, paste0(safe_name(nm), ".csv"))
      export_df <- val
      if (nrow(export_df) > 5000) {
        export_df <- export_df[seq_len(5000), , drop = FALSE]
        info$csv_truncated <- TRUE
      }
      write.csv(export_df, file = out_csv, row.names = FALSE)
      info$csv_path <- basename(out_csv)
    } else if (is.list(val) && !is.null(names(val))) {
      info$names <- names(val)
    } else if (is.vector(val)) {
      info$length <- length(val)
    }
    info
  })

  inventory <- list(
    input = input_path,
    status = "success",
    rdata_member = resolved$member,
    zip_members = resolved$zip_members,
    objects = object_info,
    api_run = NULL
  )
  write_json_simple(inventory, file.path(output_dir, "rdata_inventory.json"))

  if (resolved$cleanup && file.exists(resolved$rdata)) unlink(resolved$rdata)
  quit(status = 0)
}, error = function(e) {
  inventory <- list(
    input = input_path,
    status = "extraction_failed",
    error = conditionMessage(e),
    api_run = NULL
  )
  try(write_json_simple(inventory, file.path(output_dir, "rdata_inventory.json")), silent = TRUE)
  cat("inspect_dropviz_rdata.R error: ", conditionMessage(e), "\n", file = stderr())
  quit(status = 1)
})
