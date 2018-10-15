chomp();
@_=split(/\t/);
$acc=@_[14];
$acc =~ s/[\r\n]//g;
@traits=split(/\s*,\s*/,@_[12]);
@uris=split(/[,\s]+/,@_[13]);
for ($i=0; $i<=$#traits; ++$i)
{
  print "$acc\t$traits[$i]\t$uris[$i]\n";
}
