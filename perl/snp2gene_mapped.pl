### Input gwascat_assn.tsv, from gwascat_assn.R 
# 15. MAPPED_GENE (symbol[s])
# 16. UPSTREAM_GENE_ID (ENSG)
# 17. DOWNSTREAM_GENE_ID (ENSG)
# 18. SNP_GENE_IDS (comma delimited ENSGs)
# 22. SNPS (semicolon delimited SNPs)
# 37. STUDY_ACCESSION
###
chomp();
@_=split(/\t/);
chomp($symbs=@_[14]);
@ensgs=split(/[,\s]+/, @_[17]);
@snps=split(/[;\s]+/, @_[21]);
chomp($up_ensg=@_[15]);
chomp($down_ensg=@_[16]);
chomp($acc=@_[36]);
$acc =~ s/[\r\n]//g;
#
if ($up_ensg =~ /^ENSG/)
{
  foreach $snp (@snps)
  {
    print "$acc\t$snp\t$symbs\t$up_ensg\tmu\n";
  }
}
elsif ($down_ensg =~ /^ENSG/)
{
  foreach $snp (@snps)
  {
    print "$acc\t$snp\t$symbs\t$up_ensg\tmd\n";
  }
}
else
{
  foreach $ensg (@ensgs)
  {
    foreach $snp (@snps)
    {
      print "$acc\t$snp\t$symbs\t$ensg\tm\n";
    }
  }
}
