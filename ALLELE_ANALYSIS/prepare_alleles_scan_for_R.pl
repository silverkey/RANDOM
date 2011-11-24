#!/usr/bin/perl
use strict;
use warnings;

my $file = 'FC_alleles_analysis_CGI2/alleles.tab';
open(IN,$file);

my $start = -490;
my $pos = 0;
my $previous = 0;
my $href = {};

while(my $row = <IN>) {
  chomp($row);
  my @f = split(/\t/,$row);
  my $group = $f[0];
  my $value = $f[1];
  my @g = split(/\_/,$group);
  my $n = $g[1];
  my $t = $g[2];
  if($n == $previous) {
    print "$n\t$t\t$pos\t$value\n";
    $href->{$n}->{$pos} = $value;
    $pos += 10;
  }
  else {
    print "$n\t$t\t$start\t$value\n";
    $href->{$n}->{$start} = $value;
    $pos = $start + 10;
    $previous = $n;
  }
}
