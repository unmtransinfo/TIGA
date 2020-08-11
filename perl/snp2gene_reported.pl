### Input gwascat_assn.tsv, from gwascat_assn.R
chomp();
@_=split(/\t/);
@gsyms=split(/[,\s]+/,@_[13]);
@snps=split(/[;\s]+/,@_[21]);
$acc=@_[36];
$acc =~ s/[\r\n]//g;
foreach $gsym (@gsyms) {
  foreach $snp (@snps) {
    print "$acc\t$gsym\t$snp\tr\n";
  }
}
