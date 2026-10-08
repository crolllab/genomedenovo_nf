#!/usr/bin/env Rscript
# Plot an isoform phylogeny (IQ-TREE Newick) as a PDF.
#
# Usage:
#   plot_isoform_tree.R <tree.newick> <isoforms.fasta> <out.pdf> <seq_type> [title]
#
#   tree.newick     IQ-TREE .contree (bootstrap support as node labels) or .treefile;
#                   tips are named cds_isoform_N / protein_isoform_N
#   isoforms.fasta  ungapped isoform FASTA from summarize_isoforms.py, headers like
#                   '>cds_isoform_1 n_samples=3 status=full_length,reference'
#   seq_type        cds_isoform or protein_isoform (used in the title)
#
# The tree is midpoint-rooted. Tips read 'Isoform N (n = samples)'; isoforms
# carried only by truncated or stop-lost (extended) sequences are coloured,
# the query reference is marked with a triangle and bold label. Bootstrap
# support is shown at internal nodes when available.

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_isoform_tree.R <tree> <isoforms.fasta> <out.pdf> <seq_type> [title]")
}
tree_file <- args[1]
fasta_file <- args[2]
out_pdf <- args[3]
seq_type <- args[4]
title <- if (length(args) >= 5) args[5] else ""

col_status <- c(full_length = "black", truncated = "#D55E00", extended = "#0072B2")

# ---------------------------------------------------------------------------
# Isoform annotation from the FASTA headers
# ---------------------------------------------------------------------------
read_isoform_info <- function(path) {
  hdr <- sub("^>", "", grep("^>", readLines(path), value = TRUE))
  name <- sub("\\s.*$", "", hdr)
  field <- function(key) {
    ifelse(grepl(paste0(key, "="), hdr), sub(paste0("^.*", key, "=(\\S*).*$"), "\\1", hdr), NA)
  }
  n <- suppressWarnings(as.integer(field("n_samples")))
  statuses <- strsplit(ifelse(is.na(field("status")), "NA", field("status")), ",")
  is_ref <- vapply(statuses, function(s) "reference" %in% s, logical(1))
  # Colour class: truncated / extended when no intact sequence carries the isoform
  class <- vapply(statuses, function(s) {
    s <- setdiff(s, "reference")
    if (length(s) == 0) return("full_length")
    if (all(s %in% c("truncated", "extended"))) {
      return(if (all(s == "extended")) "extended" else "truncated")
    }
    "full_length"
  }, character(1))
  class[is_ref] <- "full_length"       # the reference is intact by definition
  data.frame(name = name, n = n, reference = is_ref, class = class,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Midpoint rooting with ape only: root on the edge holding the midpoint of the
# longest tip-to-tip path, then split that edge so both sides reach equally far.
# ---------------------------------------------------------------------------
midpoint_root <- function(tr) {
  if (Ntip(tr) < 3 || is.null(tr$edge.length)) return(tr)
  d <- cophenetic.phylo(tr)
  ij <- which(d == max(d), arr.ind = TRUE)[1, ]
  a <- ij[1]; b <- ij[2]
  half <- d[a, b] / 2
  path <- nodepath(tr, a, b)
  dist_a <- dist.nodes(tr)[a, path]
  k <- which(dist_a >= half)[1]
  if (is.na(k) || k < 2) return(tr)
  # The edge between path[k-1] and path[k]; root above its child node
  e <- which((tr$edge[, 1] == path[k - 1] & tr$edge[, 2] == path[k]) |
             (tr$edge[, 1] == path[k] & tr$edge[, 2] == path[k - 1]))
  child <- tr$edge[e, 2]
  outgroup <- if (child <= Ntip(tr)) tr$tip.label[child] else extract.clade(tr, child)$tip.label
  rt <- root(tr, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)

  # Split the combined root edges: e1 + h1 = e2 + h2 and e1 + e2 = L
  root_node <- Ntip(rt) + 1
  re <- which(rt$edge[, 1] == root_node)
  if (length(re) != 2) return(rt)
  L <- sum(rt$edge.length[re])
  depth <- node.depth.edgelength(rt)
  h <- vapply(re, function(i) {
    ch <- rt$edge[i, 2]
    tips <- if (ch <= Ntip(rt)) ch else which(rt$tip.label %in% extract.clade(rt, ch)$tip.label)
    max(depth[tips]) - depth[ch]
  }, numeric(1))
  e1 <- min(max((L + h[2] - h[1]) / 2, 0), L)
  rt$edge.length[re] <- c(e1, L - e1)
  rt
}

tr <- read.tree(tree_file)
if (is.null(tr) || Ntip(tr) < 2) stop("No tree with at least two tips in ", tree_file)
tr <- ladderize(midpoint_root(tr), right = FALSE)

info <- read_isoform_info(fasta_file)
m <- match(tr$tip.label, info$name)
num <- sub("^.*_isoform_", "", tr$tip.label)
n_lab <- ifelse(is.na(m) | is.na(info$n[m]), "?", info$n[m])
is_ref <- !is.na(m) & info$reference[m]
tip_class <- ifelse(is.na(m), "full_length", info$class[m])
labels <- paste0("Isoform ", num, " (n = ", n_lab, ifelse(is_ref, ", reference", ""), ")")
tr$tip.label <- labels

# ---------------------------------------------------------------------------
# Plot: page height scales with the number of tips
# ---------------------------------------------------------------------------
n_tips <- Ntip(tr)
max_depth <- max(node.depth.edgelength(tr))
pdf(out_pdf, width = 8, height = max(4, 0.25 * n_tips + 2.5))
par(mar = c(1, 1, 3, 1))

# Extra space below the tips holds the scale bar and the legend
plot.phylo(tr, type = "phylogram", show.tip.label = TRUE, cex = 0.8,
           tip.color = col_status[tip_class], font = ifelse(is_ref, 2, 1),
           label.offset = 0.02 * max_depth, no.margin = FALSE,
           x.lim = c(0, max_depth * 1.7), y.lim = c(-2, n_tips + 0.5))

if (any(is_ref)) {
  tiplabels(tip = which(is_ref), pch = 17, col = "black", cex = 0.9, offset = 0)
}

# Bootstrap support on internal nodes (IQ-TREE writes it as node labels)
if (!is.null(tr$node.label)) {
  support <- suppressWarnings(as.numeric(sub("/.*$", "", tr$node.label)))
  show <- setdiff(which(!is.na(support)), 1)    # no support value for the (midpoint) root
  if (length(show) > 0) {
    nodelabels(text = round(support[show]), node = Ntip(tr) + show, frame = "none",
               adj = c(1.2, -0.4), cex = 0.6, col = "grey30")
  }
}

add.scale.bar(x = 0, y = -0.2, cex = 0.7)
type_label <- if (grepl("^protein", seq_type)) "Protein isoforms" else "CDS isoforms"
main <- if (nzchar(title)) paste0(type_label, ": ", title) else type_label
title(main = main, cex.main = 0.9, font.main = 1)

legend(x = 0, y = -0.9, yjust = 0.5, horiz = TRUE, bty = "n", cex = 0.7,
       legend = c("full-length", "truncated only", "stop codon lost only", "query reference"),
       text.col = c(col_status[c("full_length", "truncated", "extended")], "black"),
       pch = c(NA, NA, NA, 17), col = "black")
invisible(dev.off())
