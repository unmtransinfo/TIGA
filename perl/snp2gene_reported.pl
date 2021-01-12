### Input gwascat_assn.tsv, from gwascat_assn.R
chomp();
@_=split(/\t/);
$acc=@_[36];
$acc =~ s/[\r\n]//g;
@snps=split(/[;\s]+/,@_[21]);
@gsyms=split(/[,\s]+/,@_[13]);
foreach $gsym (@gsyms) {
  foreach $snp (@snps) {
    print "$acc\t$snp\t$gsym\t\tr\n";
  }
}
