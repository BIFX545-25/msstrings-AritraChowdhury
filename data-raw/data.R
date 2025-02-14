#create sample data set for package
library(UniProt.ws)
library(Biostrings)
library(stringr)

seqs <- queryUniProt('P51681', 'sequence')$Sequence

AAStringSet(seqs)

#If we are using Trypsin in our Mass Spec protocol, we’ll expect this to
#be split at either lysine (K) or arginine (R), except when followed by a
#proline (P). So, these are the peptide sequences we might expect to get
#in our MS results:

# cut points (plus the ends of the sequence)
cp <- c(0, str_locate_all(seqs, '[R|K]')[[1]][,'start'], nchar(seqs))

# possible fragments
frags <- str_sub(seqs, cp[-length(cp)] + 1, cp[-1])

# Glycosylation
frags[1] <- str_replace(frags[1], 'S', 'S[Glycosylation (S)]')

# Oxidation
frags[2] <- str_replace(frags[2], 'M', 'M[Oxidation (M)]')

# Acetylation
frags[3] <- str_replace(frags[3], 'Y', 'Y[Acetylation (Y)]')

# add a second modification (Deamidation) to #3
frags[3] <- str_replace(frags[3], 'Q', 'Q[Deamidation (Q)]')

# add underscores as Spectronaut does
frags <- paste0('_', frags, '_')

MSfrags<- frags
usethis::use_data(MSfrags, overwrite=TRUE)