#!/usr/bin/env Rscript

library(readr)
library(data.table)

message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
ODIR <- "data/20210212_GCAPI"

ifile_file <- ifelse(length(args)>0, args[1], paste0(ODIR, "/gwascat_snp2gene_FILE.tsv"))
ifile_api <- ifelse(length(args)>1, args[2], paste0(ODIR, "/gwascat_snp2gene_API.tsv"))
ofile	<- ifelse(length(args)>2, args[3], paste0(ODIR, "/gwascat_snp2gene_MERGED.tsv"))

message(sprintf("Input FILE file: %s", ifile_file))
message(sprintf("Input API file: %s", ifile_api))
message(sprintf("Output file: %s", ofile))

snp2gene_file <- read_delim(ifile_file, "\t", col_types=cols(.default=col_character()))
setDT(snp2gene_file )
message(sprintf("FILE (%s) rows: %d", sub("^.*/", "", ifile_file), nrow(snp2gene_file)))
snp2gene_file <- snp2gene_file[MAPPED_OR_REPORTED != "r"] #mapped only (m, mu, md)
setnames(snp2gene_file, old=c("SNPS"), new=c("SNP"))
message(sprintf("FILE unique SNPs: %d", snp2gene_file[, uniqueN(SNP)]))
message(sprintf("FILE unique ENSGs: %d", snp2gene_file[, uniqueN(ENSG)]))
message(sprintf("FILE unique SNP2GENE pairs: %d", nrow(unique(snp2gene_file[, .(SNP, ENSG)]))))

snp2gene_api <- read_delim(ifile_api, "\t", col_types=cols(.default=col_character(), isIntergenic=col_logical(), isUpstream=col_logical(), isDownstream=col_logical(), isClosestGene=col_logical(), distance=col_integer()))
setDT(snp2gene_api )
message(sprintf("API (%s) rows: %d", sub("^.*/", "", ifile_api), nrow(snp2gene_api)))
setnames(snp2gene_api, old=c("rsId", "geneName"), new=c("SNP", "GSYMB"))
snp2gene_api[, ENSG := ensemblGeneIds]
snp2gene_api[, MAPPED_OR_REPORTED := ifelse((isUpstream), "mu", ifelse((isDownstream), "md", "m"))]
snp2gene_api <- snp2gene_api[!is.na(ENSG)]
snp2gene_api <- unique(snp2gene_api[, .(SNP, GSYMB, ensemblGeneIds, ENSG, MAPPED_OR_REPORTED)])
setkey(snp2gene_api, "SNP")
snp2gene_api <- unique(snp2gene_api[, list(ENSG=unlist(strsplit(ENSG, "[ ,]"))), by=c("SNP", "GSYMB", "ensemblGeneIds", "MAPPED_OR_REPORTED")])
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

#Check for our favorite gene (SLC25A44), SNPs, studies:
snp2gene_test <- snp2gene[ENSG=="ENSG00000160785" & SNP %in% c("rs2273833", "rs6684514", "rs144991356")]
message(sprintf("Check for SLC25A44 (ENSG00000160785): %d rows", nrow(snp2gene[ENSG=="ENSG00000160785"])))
gcsts_test <- c("GCST002390", "GCST005145", "GCST006001")
message(sprintf("Check HbA1c measurement (EFO_0004541) GCSTs (%s): %d rows", paste(collapse=",", gcsts_test), nrow(snp2gene[STUDY_ACCESSION %in% gcsts_test])))
snps_test <- c("rs2273833", "rs6684514", "rs144991356")
message(sprintf("Check HbA1c measurement (EFO_0004541) SNPS (%s): %d rows", paste(collapse=",", snps_test), nrow(snp2gene[SNP %in% snps_test])))
message(sprintf("Check for SLC25A44 - HbA1c measurement (ENSG00000160785-EFO_0004541): %s", nrow(snp2gene_test)>0))
print(snp2gene_test)
#
write_delim(snp2gene, ofile, delim="\t")
