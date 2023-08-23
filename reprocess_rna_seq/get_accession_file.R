project <- "project_name" # Example: PRJNA315611 for GSE79362
url <- sprintf("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=%s&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true&limit=0", 
               project)
destfile <- "path_to_file" # Example: "~/Downloads/prjna.txt"
download.file(url, destfile)
prjna <- read.delim(destfile)
geo_accession <- prjna$run_accession
write.table(prjna$run_accession, "geo_accession.txt",
            quote = F, row.names = F, col.names = F)

