#Description: Create a blast database for the protein secretory pathway 
#Author: Lisa Yeter Mesrop 
#Goal: Use the blast database to find protein secretory pathway genes in the significantly upregulated genes of the bioluminescent upper lip, non-luminous upper lip, and BCN. 

#download the protein secretory pathway genes from Feizi et al. (2017)
Hu_secretory_system_genes_ensembleids <- as.character(Hu_secretory_system_genes$ensgid)

#rename
ensembl_ids <- Hu_secretory_system_genes_ensembleids

# function to fetch sequence from Ensembl REST API
fetch_sequence <- function(ensembl_id) {
  url <- paste0("https://rest.ensembl.org/sequence/id/", ensembl_id, "?content-type=text/x-fasta")
  response <- GET(url)
  
  if (status_code(response) != 200) {
    warning("Failed to fetch sequence for ID: ", ensembl_id)
    return(NULL)
  }
  
  content(response, "text")
}

# Fetch sequences and filter them
sequences <- sapply(ensembl_ids, fetch_sequence)

# Remove NULL values (in case of failed fetches)
sequences <- sequences[!sapply(sequences, is.null)]

# Combine sequences into a single character vector
combined_sequences <- paste(sequences, collapse = "\n")

# Write sequences to output file
output_file <- "secretory_protein_pathway_fasta.txt"
writeLines(combined_sequences, con = output_file)


#make blast datatbase
makeblastdb -in secretory_protein_pathway.fasta -dbtype nucl -out secretory_protein_pathway_ensembl_db

#references
Feizi, A., Gatto, F., Uhlen, M. et al. Human protein secretory pathway genes are expressed in a tissue-specific pattern to match processing demands of the secretome. npj Syst Biol Appl 3, 22 (2017)
