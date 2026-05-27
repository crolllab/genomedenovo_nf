library(jsonlite)

# Keys from summary.before_filtering / summary.after_filtering
SUMMARY_KEYS <- c("total_reads", "total_bases", "q20_bases", "q30_bases",
                  "q20_rate", "q30_rate", "read1_mean_length",
                  "read2_mean_length", "gc_content")

# Keys from filtering_result (appear once)
FILTER_KEYS <- c("passed_filter_reads", "low_quality_reads",
                 "too_many_N_reads", "too_short_reads", "too_long_reads")

files <- list.files(pattern = "_(PE|SE)_fastp.json$")
if (length(files) == 0) quit(save = "no", status = 0)

samples <- sub("_(PE|SE)_fastp.json$", "", files)

# Row layout: before_filtering keys, then after_filtering keys (same names
# prefixed for clarity), then filtering_result keys, then extra metrics.
before_row_names <- paste0("before_", SUMMARY_KEYS)
after_row_names  <- paste0("after_",  SUMMARY_KEYS)
row_names <- c(before_row_names, after_row_names, FILTER_KEYS,
               "duplication_rate", "insert_size_peak")

mat <- matrix(NA_real_, nrow = length(row_names), ncol = length(files),
              dimnames = list(row_names, samples))

n_sum    <- length(SUMMARY_KEYS)
n_filter <- length(FILTER_KEYS)

for (i in seq_along(files)) {
  j  <- fromJSON(files[[i]], simplifyVector = TRUE)
  bf <- j$summary$before_filtering
  af <- j$summary$after_filtering
  fr <- j$filtering_result

  for (k in seq_along(SUMMARY_KEYS)) {
    key <- SUMMARY_KEYS[[k]]
    mat[k,          i] <- if (!is.null(bf[[key]])) as.numeric(bf[[key]]) else NA_real_
    mat[k + n_sum,  i] <- if (!is.null(af[[key]])) as.numeric(af[[key]]) else NA_real_
  }
  for (k in seq_along(FILTER_KEYS)) {
    key <- FILTER_KEYS[[k]]
    mat[2*n_sum + k, i] <- if (!is.null(fr[[key]])) as.numeric(fr[[key]]) else NA_real_
  }
  base <- 2*n_sum + n_filter
  mat[base + 1, i] <- tryCatch({
    v <- as.numeric(j$duplication$rate)
    if (length(v) == 0) NA_real_ else v
  }, error = function(e) NA_real_)
  mat[base + 2, i] <- tryCatch({
    v <- as.numeric(j$insert_size$peak)
    if (length(v) == 0) NA_real_ else v
  }, error = function(e) NA_real_)
}

# Write TSV (metrics as rows, samples as columns)
write.table(mat, file = "fastp_summary.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Write human-readable TXT
con <- file("fastp_summary.txt", open = "wt")

header_line <- paste(rep("=", 65), collapse = "")
writeLines(header_line, con)
writeLines("  FASTP Trimming Statistics", con)
writeLines(paste0("  Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), con)
writeLines(header_line, con)
writeLines("", con)

# Determine column widths
metric_width <- max(nchar(row_names)) + 2
val_width    <- max(10, max(nchar(samples)) + 2)

# Header row
header <- formatC("Metric", width = -metric_width, flag = "-")
for (s in samples) header <- paste0(header, formatC(s, width = val_width, flag = " "))
writeLines(header, con)
writeLines(paste(rep("-", nchar(header)), collapse = ""), con)

# Data rows
for (rn in row_names) {
  row_str <- formatC(rn, width = -metric_width, flag = "-")
  for (s in samples) {
    val <- mat[rn, s]
    if (is.na(val)) {
      cell <- "NA"
    } else if (rn %in% c("before_q20_rate", "before_q30_rate", "before_gc_content",
                          "after_q20_rate",  "after_q30_rate",  "after_gc_content",
                          "duplication_rate")) {
      cell <- sprintf("%.4f", val)
    } else {
      cell <- format(round(val), big.mark = ",", scientific = FALSE)
    }
    row_str <- paste0(row_str, formatC(cell, width = val_width, flag = " "))
  }
  writeLines(row_str, con)
}

writeLines("", con)
writeLines(header_line, con)
close(con)
