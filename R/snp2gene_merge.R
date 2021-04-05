#!/usr/bin/env Rscript

library(readr)
library(data.table)

message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
ODIR <- "data/20210212_GCAPI"

ifile_file <- ifelse(length(args)>0, args[1], paste0(ODIR, "/gwascat_snp2gene_FILE.tsv"))
ifile_api <- ifelse(length(args)>1, args[2], paste0(ODIR, "/gwascat_snp2gene_API.tsv"))
ofile	<- ifelse(length(args)>2, args[3], paste0(ODIR, "/gwascat_snp2gene_MERGED.tsv"))

if (length(args)>3) {
  message("ERROR: Syntax: snp2gene_merge.R [SNP2GENE_FILE_FILE [SNP2GENE_API_FILE [OFILE]]]\n...or... no args for defaults")
  quit()
}

writeLines(sprintf("Input FILE file: %s", ifile_file))
writeLines(sprintf("Input API file: %s", ifile_api))
writeLines(sprintf("Output file: %s", ofile))

snp2gene_file <- read_delim(ifile_file, "\t", col_types=cols(.default=col_character()))
setDT(snp2gene_file )
message(sprintf("FILE (%s) rows: %d", sub("^.*/", "", ifile_file), nrow(snp2gene_file)))
snp2gene_file <- snp2gene_file[REPORTED_OR_MAPPED != "r"] #mapped only (m, mu, md)
message(sprintf("FILE unique SNPs: %d", snp2gene_file[, uniqueN(SNP)]))
message(sprintf("FILE unique ENSGs: %d", snp2gene_file[, uniqueN(ENSG)]))
message(sprintf("FILE unique SNP2GENE pairs: %d", nrow(unique(snp2gene_file[, .(SNP, ENSG)]))))

snp2gene_api <- read_delim(ifile_api, "\t", col_types=cols(.default=col_character(), isIntergenic=col_logical(), isUpstream=col_logical(), isDownstream=col_logical(), isClosestGene=col_logical(), distance=col_integer()))
setDT(snp2gene_api )
message(sprintf("API (%s) rows: %d", sub("^.*/", "", ifile_api), nrow(snp2gene_api)))
setnames(snp2gene_api, old=c("rsId", "geneName"), new=c("SNP", "GSYMB"))
snp2gene_api[, ENSG := ensemblGeneIds]
snp2gene_api[, REPORTED_OR_MAPPED := ifelse((isUpstream), "mu", ifelse((isDownstream), "md", "m"))]
snp2gene_api <- snp2gene_api[!is.na(ENSG)]
snp2gene_api <- unique(snp2gene_api[, .(SNP, GSYMB, ensemblGeneIds, ENSG, REPORTED_OR_MAPPED)])
setkey(snp2gene_api, "SNP")
snp2gene_api <- unique(snp2gene_api[, list(ENSG=unlist(strsplit(ENSG, "[ ,]"))), by=c("SNP", "GSYMB", "ensemblGeneIds", "REPORTED_OR_MAPPED")])
snp2gene_api[, ensemblGeneIds := NULL]
#Get STUDY_ACCESSIONs from FILE file.
snp2gene_api <- unique(merge(snp2gene_api, snp2gene_file[, .(STUDY_ACCESSION, SNP)], by="SNP", allow.cartesian=T))
message(sprintf("API unique SNPs: %d", snp2gene_api[, uniqueN(SNP)]))
message(sprintf("API unique ENSGs: %d", snp2gene_api[, uniqueN(ENSG)]))
message(sprintf("API unique SNP2GENE pairs: %d", nrow(unique(snp2gene_api[, .(SNP, ENSG)]))))

message(sprintf("API unique SNPs not in FILE: %d",
	length(setdiff(snp2gene_api[, unique(SNP)], snp2gene_file[, unique(SNP)]))))
message(sprintf("API unique ENSGs not in FILE: %d",
	length(setdiff(snp2gene_api[, unique(ENSG)], snp2gene_file[, unique(ENSG)]))))
message(sprintf("FILE unique SNPs not in API: %d",
	length(setdiff(snp2gene_file[, unique(SNP)], snp2gene_api[, unique(SNP)]))))
message(sprintf("FILE unique rs-SNPs not in API: %d",
                sum(grepl("^rs", setdiff(snp2gene_file[, unique(SNP)], snp2gene_api[, unique(SNP)])))))
message(sprintf("FILE unique ENSGs not in API: %d",
	length(setdiff(snp2gene_file[, unique(ENSG)], snp2gene_api[, unique(ENSG)]))))

write_delim(data.table(SNP=setdiff(snp2gene_file[, unique(SNP)], snp2gene_api[, unique(SNP)])), paste0(ODIR, "/gwascat_snp2gene_FILE_NOT_API.rs"), delim="\t")

snp2gene <- rbindlist(list(snp2gene_file, snp2gene_api), use.names=T)
#snp2gene <- snp2gene[grepl("^rs", SNP)] # Should we keep non-rs SNPs?
message(sprintf("TOTAL SNP2GENE pairs: %d", nrow(unique(snp2gene[, .(SNP, ENSG)]))))

#Check for our favorite gene (SLC25A44) and SNPs:
snp2gene_test <- snp2gene[ENSG=="ENSG00000160785" & SNP %in% c("rs2273833", "rs6684514", "rs144991356")]
message(sprintf("Check for SLC25A44 - HbA1c measurement (ENSG00000160785-EFO_0004541): %s", nrow(snp2gene_test)>0))
print(snp2gene_test)

write_delim(snp2gene, ofile, delim="\t")
