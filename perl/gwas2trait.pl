# Problem: Commas may exist in single trait names.
chomp();
@_=split(/\t/);
$acc=@_[14];
$acc =~ s/[\r\n]//g;
@traits=split(/\s*,\s*/, @_[12]);
@uris=split(/\s*,\s*/, @_[13]);
if ($#uris > 1) {
  if ($#traits == $#uris) {
    for ($i=0; $i<=$#uris; ++$i) {
      print "$acc\t$traits[$i]\t$uris[$i]\n";
    }
  }
  else {
    for ($i=0; $i<=$#uris; ++$i) {
      print "$acc\t@_[12]\t$uris[$i]\n";
    }
  }
}
else {
  print "$acc\t@_[12]\t@_[13]\n";
}
