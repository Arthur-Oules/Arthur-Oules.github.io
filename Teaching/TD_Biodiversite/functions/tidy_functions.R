Get_sequence <- function(AC_number) {
# Gets nucleotide sequences from Genbank Accession number
  if (is.na(AC_number) == TRUE) {NA} else {
    entrez_fetch(
      db      = "nuccore",
      id      = AC_number,
      rettype = "fasta_cds_na",
      retmode = "text"
    ) |>
      strsplit("]") |>
      unlist() |>
      tail(1) |>
      gsub(pattern = "\n", replacement = "",  x = _)
  }
}

write_to_fasta <- function(tib, path) {
# Writes two column tib as a .fa file
# Tibble structure must be :
# fasta name | sequence
  tib$accession <- paste0(">", tib$accession)
  c(rbind(tib$accession, tib$sequences)) |> write(file = path)
}