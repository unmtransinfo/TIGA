### Input gwascat_assn.tsv, from gwascat_assn.R and gwas_catalog_v1.0.2-associations_e94_r2018-09-30.tsv
chomp();
@_=split(/\t/);
@gsyms=split(/[;,]/,@_[14]);
@snps=split(/[;\s]+/,@_[21]);
chomp($acc=@_[36]);
$acc =~ s/[\r\n]//g;
foreach $gsym (@gsyms)
{
  $gsym =~ s/^ *(\S+)\s*/$1/g;
  if ($gsym =~ /-/)
  {
    ($u,$d)=split(/ *- */,$gsym);
    $u =~ s/^ *(\S+)\s*/$1/g;
    $d =~ s/^ *(\S+)\s*/$1/g;
    foreach $snp (@snps)
    {
      print "$acc\t$u\t$snp\tmu\n";
      print "$acc\t$d\t$snp\tmd\n";
    }
  }
  elsif ($gsym =~ / x /)
  {
    ($u,$d)=split(/ +x +/,$gsym);
    $u =~ s/^ *(\S+)\s*/$1/g;
    $d =~ s/^ *(\S+)\s*/$1/g;
    foreach $snp (@snps)
    {
      print "$acc\t$u\t$snp\tmu\n";
      print "$acc\t$d\t$snp\tmd\n";
    }
  }
  else
  {
    foreach $snp (@snps)
    {
      print "$acc\t$gsym\t$snp\tm\n";
    }
  }
}
