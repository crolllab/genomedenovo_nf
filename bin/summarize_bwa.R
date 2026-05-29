library(jsonlite)

# Keys to extract from the "QC-passed reads" section of samtools flagstat JSON
KEYS <- c(
  "total", "primary", "secondary", "supplementary",
  "mapped", "mapped %", "primary mapped", "primary mapped %",
  "paired in sequencing", "properly paired", "properly paired %",
  "singletons", "singletons %"
)

PERCENT_KEYS <- c("mapped %", "primary mapped %", "properly paired %", "singletons %")

files <- list.files(pattern = "_flagstat\\.json$")
if (length(files) == 0) quit(save = "no", status = 0)

samples <- sub("_flagstat\\.json$", "", files)

mat <- matrix(NA_real_, nrow = length(KEYS), ncol = length(files),
              dimnames = list(KEYS, samples))

for (i in seq_along(files)) {
  j   <- fromJSON(files[[i]], simplifyVector = TRUE)
  qc  <- j[["QC-passed reads"]]
  for (k in seq_along(KEYS)) {
    val <- qc[[KEYS[[k]]]]
    mat[k, i] <- if (!is.null(val) && length(val) > 0 && !is.na(val)) as.numeric(val) else NA_real_
  }
}

# Write TSV (metrics as rows, samples as columns)
write.table(mat, file = "bwa_mapping_summary.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Write human-readable TXT
con <- file("bwa_mapping_summary.txt", open = "wt")

header_line <- paste(rep("=", 65), collapse = "")
writeLines(header_line, con)
writeLines("  BWA Mapping Statistics (QC-passed reads)", con)
writeLines(paste0("  Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), con)
writeLines(header_line, con)
writeLines("", con)

# Determine column widths
metric_width <- max(nchar(KEYS)) + 2
val_width    <- max(10, max(nchar(samples)) + 2)

# Header row
header <- formatC("Metric", width = -metric_width, flag = "-")
for (s in samples) header <- paste0(header, formatC(s, width = val_width, flag = " "))
writeLines(header, con)
writeLines(paste(rep("-", nchar(header)), collapse = ""), con)

# Data rows
for (rn in KEYS) {
  row_str <- formatC(rn, width = -metric_width, flag = "-")
  for (s in samples) {
    val <- mat[rn, s]
    if (is.na(val)) {
      cell <- "NA"
    } else if (rn %in% PERCENT_KEYS) {
      cell <- sprintf("%.2f%%", val)
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
