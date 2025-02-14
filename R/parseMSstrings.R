# parseMSstrings.R

#' parse_mods
#' Utilities for parsing MS strings
#' 
#' This function takes strings from Spectronaut output, pulls metadata on protein modifications, and displays
#' this information in a nice format
#' 
#' @param seqs Character vector, peptide sequence strings from Spectronaut output
#' @param format Character value, denotes the format of input strings
#' 
#' @return A `protein_mods` object, which is a data fram with the animo acid sequences and modifications.
#' @export
#' @importFrom Biostrings AAStringSet
#' @importFrom dplyr tibble
#' @importFrom purrr map map_chr
#' @importFrom stringr fixed str_extract_all str_replace str_replace_all str_split


# library(Biostrings) # for `AAStringSet` - BiocManager::install('Biostrings')
# library(dplyr) # for `tibble`
# library(purrr) # for `map` and `map_chr`
# library(stringr) # for `str_extract_all`, `str_replace`, `str_replace_all`, and `str_split`

parse_mods <- function(seqs, format = 'Spectronaut')
{
  # remove '_'
  seqs <- seqs |>
    str_replace_all('_', '')


  # locate modification(s)
  mods <- str_extract_all(seqs, "\\[.*?\\]") |>
    map(~ .x |>
          str_replace(fixed("["), '') |>
          str_replace(fixed("]"), '') |>
          str_replace("\\(.*?\\)", '') |>
          trimws())


  # remove modifications from seqs
  seqs <- str_split(seqs, "\\[.*?\\]")


  # find where the mods are at
  mods_at <- map(seqs, ~
                   {
                     retval <- nchar(.x) |>
                       cumsum()

                     retval[-length(retval)] # this is the end of the peptide, not a modification
                   })


  # return the data.frame
  retval <- tibble(sequence = map_chr(seqs, ~ paste(.x, collapse = '')) |>
                     AAStringSet() |>
                     as.list(),
                   mods = mods,
                   mods_at = mods_at)

  class(retval) <- c('protein_mods', class(retval))

  return(retval)
}
