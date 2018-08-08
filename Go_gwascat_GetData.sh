#!/bin/sh
#############################################################################
### Go_gwascat_GetData.sh - Create CSV files for import to RDB:
###
### GWAS Catalog from http://www.genome.gov/gwastudies/
### Migrated to: http://www.ebi.ac.uk/gwas/
#############################################################################
### GWAS Catalog studies each have a "study_accession".
### Also are associated with a publication (PubMedID), but not uniquely.
### 
### See https://www.ebi.ac.uk/gwas/docs/fileheaders
### 
### MAPPED GENE(S): Gene(s) mapped to the strongest SNP. If the SNP is located
### within a gene, that gene is listed. If the SNP is intergenic, the upstream
### and downstream genes are listed, separated by a hyphen.
#############################################################################
### "OR or BETA: Reported odds ratio or beta-coefficient associated with
### strongest SNP risk allele. Note that if an OR <1 is reported this is
### inverted, along with the reported allele, so that all ORs included in
### the Catalog are >1. Appropriate unit and increase/decrease are included
### for beta coefficients."
#############################################################################
### Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H,
### Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H. The NHGRI
### GWAS Catalog, a curated resource of SNP-trait associations. Nucleic
### Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006.
#############################################################################
### iCite annotations from iCite API, with all PMIDs from GWASCatalog.
#############################################################################
### Jeremy Yang
###  3 Nov 2017
#############################################################################
#
SRCDATADIR="/home/data/gwascatalog/data"
DATADIR="data"
#
tempperl="perl/tmp.pl"
#
#Source files:
#gwasfile="${SRCDATADIR}/gwas_catalog_v1.0.1-studies_r2017-02-13.tsv"
gwasfile="${SRCDATADIR}/gwas_catalog_v1.0.1-studies_r2017-10-10.tsv"
#
#assnfile="${SRCDATADIR}/gwas_catalog_v1.0.1-associations_e87_r2017-02-13.tsv"
assnfile="${SRCDATADIR}/gwas_catalog_v1.0.1-associations_e90_r2017-10-10.tsv"
#
#Cleaned output files:
csvfile_assn="${DATADIR}/gwascat_assn.csv"
#
#Clean non-ascii (Greek)
tsvfile_gwas="${DATADIR}/gwascat_gwas.tsv"
unicode_utils.py --i $gwasfile --o $tsvfile_gwas
csv_utils.py --tsv --i $tsvfile_gwas --fixtags --overwrite_input_file
#
#Clean, separate OR_or_beta into oddsratio, beta columns:
R/gwascat_assn.R $assnfile $csvfile_assn
#
#############################################################################
### Parse initial_sample_size
#	cat <<__EOF__ >${tempperl}
#	chomp();
#	@_=split(/\t/);
#	\$acc=@_[0];
#	\$acc =~ s/[\r\n]//g;
#	\$iss = @_[1];
#	\$iss =~ s/[^ 0-9]//g;
#	\$iss =~ s/^ *//g;
#	\$iss =~ s/ *$//g;
#	\$iss =~ s/ +/+/g;
#	\$iss = eval \$iss;
#	print "\"\$acc\"\t\$iss\n";
#	__EOF__
#	#
#	csv_utils.py --tsv --i $tsvfile_gwas --subsetcols \
#		--coltags "STUDY_ACCESSION,INITIAL_SAMPLE_SIZE" \
#		|sed -e '1d' \
#		|perl -n ${tempperl}
#
#############################################################################
#trait links:
traitfile="${DATADIR}/gwascat_trait.tsv"
#
#SNP to gene links:
snp2genefile="${DATADIR}/gwascat_snp2gene.tsv"
#
#############################################################################
### TRAITS:
####
cat <<__EOF__ >${tempperl}
chomp();
@_=split(/\t/);
\$acc=@_[0];
\$acc =~ s/[\r\n]//g;
@ts=split(/\s*,\s*/,@_[1]);
@uris=split(/[,\s]+/,@_[2]);
for (\$i=0; \$i<=\$#ts; ++\$i)
{
  print "\"\$acc\"\t\"\$ts[\$i]\"\t\"\$uris[\$i]\"\n";
}
__EOF__
#
printf "study_accession\ttrait\ttrait_uri\n" >${traitfile}
#
csv_utils.py --tsv --i $tsvfile_gwas --subsetcols \
	--coltags "STUDY_ACCESSION,MAPPED_TRAIT,MAPPED_TRAIT_URI" \
	|sed -e '1d' \
	|perl -n ${tempperl} \
	>>${traitfile}
#
#############################################################################
### REPORTED GENES:
cat <<__EOF__ >${tempperl}
chomp();
@_=split(/\t/);
@gsyms=split(/[,\s]+/,@_[0]);
@snps=split(/[;\s]+/,@_[1]);
\$acc=@_[2];
\$acc =~ s/[\r\n]//g;
foreach \$gsym (@gsyms) {
  foreach \$snp (@snps) {
    print "\"\$acc\"\t\"\$gsym\"\t\"\$snp\"\tr\n";
  }
}
__EOF__
#
printf "study_accession\tgsymb\tsnp\treported_or_mapped\n" >${snp2genefile}
#
csv_utils.py --i $csvfile_assn --subsetcols \
	--coltags "REPORTED_GENE(S),SNPS,STUDY_ACCESSION" \
	|sed -e '1d' \
	|perl -n ${tempperl} \
	>>${snp2genefile}
#
#############################################################################
### MAPPED GENES:
### Separate mapped into up-/down-stream.
# "m" - mapped within gene
# "mu" - mapped to upstream gene
# "md" - mapped to downstream gene
###
#
cat <<__EOF__ >${tempperl}
chomp();
@_=split(/\t/);
@gsyms=split(/[;,]/,@_[0]);
@snps=split(/[;\s]+/,@_[1]);
chomp(\$acc=@_[2]);
\$acc =~ s/[\r\n]//g;
foreach \$gsym (@gsyms)
{
  \$gsym =~ s/^ *(\\S+)\s*/\$1/g;
  if (\$gsym =~ /-/)
  {
    (\$u,\$d)=split(/ *- */,\$gsym);
    \$u =~ s/^ *(\\S+)\s*/\$1/g;
    \$d =~ s/^ *(\\S+)\s*/\$1/g;
    foreach \$snp (@snps)
    {
      print "\"\$acc\"\t\"\$u\"\t\"\$snp\"\tmu\n";
      print "\"\$acc\"\t\"\$d\"\t\"\$snp\"\tmd\n";
    }
  }
  elsif (\$gsym =~ / x /)
  {
    (\$u,\$d)=split(/ +x +/,\$gsym);
    \$u =~ s/^ *(\\S+)\s*/\$1/g;
    \$d =~ s/^ *(\\S+)\s*/\$1/g;
    foreach \$snp (@snps)
    {
      print "\"\$acc\"\t\"\$u\"\t\"\$snp\"\tmu\n";
      print "\"\$acc\"\t\"\$d\"\t\"\$snp\"\tmd\n";
    }
  }
  else
  {
    foreach \$snp (@snps)
    {
      print "\"\$acc\"\t\"\$gsym\"\t\"\$snp\"\tm\n";
    }
  }
}
__EOF__
#
csv_utils.py --i $csvfile_assn --subsetcols \
	--coltags "MAPPED_GENE,SNPS,STUDY_ACCESSION" \
	|sed -e '1d' \
	| perl -n ${tempperl} \
	>>${snp2genefile}
#
#############################################################################
###
csv_utils.py --i data/gwascat.tsv --tsv --coltag PUBMEDID --extractcol \
	|sort -nu >data/gwascat.pmid
#
pubmed_icite.py list \
	--i data/gwascat.pmid \
	--o data/gwascat_icite.csv
#
unicode_utils.py --i data/gwascat_icite.csv --o data/gwascat_icite_ascii.csv
mv data/gwascat_icite_ascii.csv data/gwascat_icite.csv
#
#############################################################################
# gt_stats.csv generated by this R script:
#
./R/gwascat_gt_stats.R
#
