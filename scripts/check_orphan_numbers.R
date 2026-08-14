#!/usr/bin/env Rscript
# Orphan number detector: ensures all numeric values in .md files appear in NUMBERS_TABLE.tsv
# Ignores: dates, GSE accessions, DOIs, version strings, page/line numbers, ORCIDs, integers <100

library(stringr)

# Find all .md files in project root and subdirectories (exclude inst/, data/, .git/)
md_files <- list.files(
  path = ".",
  pattern = "\\.md$",
  recursive = TRUE,
  ignore.case = TRUE,
  full.names = TRUE
)
md_files <- md_files[!grepl("^\\./(inst/|data/|\\.git|vignettes/)", md_files)]

# Load NUMBERS_TABLE.tsv as character strings to preserve precision
if (!file.exists("NUMBERS_TABLE.tsv")) {
  cat("ERROR: NUMBERS_TABLE.tsv not found in project root\n")
  quit(status = 1)
}
numbers_table <- read.delim("NUMBERS_TABLE.tsv", stringsAsFactors = FALSE, header = TRUE, colClasses = "character")
valid_numbers <- unique(numbers_table[[1]])

# Extract all numeric values from .md files
orphans <- list()

for (md_file in md_files) {
  lines <- readLines(md_file, warn = FALSE)

  for (line_num in seq_along(lines)) {
    line <- lines[line_num]

    # Extract all numeric patterns
    # Pattern: numbers with optional leading minus, decimal point, scientific notation
    matches <- str_match_all(
      line,
      "(?<![\\w.-])([0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)(?![\\w.-])"
    )[[1]]

    if (nrow(matches) > 0) {
      for (i in seq_len(nrow(matches))) {
        value <- matches[i, 2]

        # Skip: dates (4-digit year patterns like 2026), GSE accessions (GSE\d+),
        # DOIs (10.XXXX), version strings (v0.1.0), page/line numbers preceded by p. or line,
        # ORCIDs (0000-XXXX-XXXX-XXXX), integers < 100

        # Skip if part of date
        if (grepl("^(19|20|21)[0-9]{2}$", value)) next

        # Skip if part of DOI or accession
        if (grepl("(GSE|GEO|PMID|Zenodo|DOI|10\\.[0-9]+)", line)) next

        # Skip if version string (e.g., v0.1.0)
        if (grepl("v[0-9]+\\.[0-9]+", line)) next

        # Skip if page/line number (preceded by p., pp., line, lines)
        if (grepl("\\b(p\\.|pp\\.|line|lines)\\s*[0-9]", line)) next

        # Skip if ORCID pattern
        if (grepl("0000-[0-9]{4}-[0-9]{4}-[0-9]{4}", line)) next

        # Skip integers < 100
        if (grepl("^[0-9]+$", value) && as.numeric(value) < 100) next

        # Check if number is in NUMBERS_TABLE
        if (!(value %in% valid_numbers)) {
          orphans[[length(orphans) + 1]] <- data.frame(
            file = md_file,
            line_num = line_num,
            value = value,
            context = substr(line, 1, 80),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

# Report results
if (length(orphans) > 0) {
  orphan_df <- do.call(rbind, orphans)
  cat("ORPHAN NUMBERS FOUND:\n")
  cat(sprintf("%.0f orphan value(s) in markdown files\n\n", nrow(orphan_df)))

  for (i in seq_len(nrow(orphan_df))) {
    cat(sprintf(
      "File: %s\nLine %d: %s\nOrphan: %s\n\n",
      orphan_df$file[i],
      orphan_df$line_num[i],
      orphan_df$context[i],
      orphan_df$value[i]
    ))
  }

  quit(status = 1)
} else {
  cat("OK: All numeric values in markdown files appear in NUMBERS_TABLE.tsv\n")
  quit(status = 0)
}
