#!/usr/bin/perl
use strict;
use warnings;

my $input = $ARGV[0];
open(IN,$input);
my $header = <IN>;
my $output = "simple_$input";
open(OUT,">$output");
print OUT join("\t",qw(id chr start end strand length tr.name tr.strand))."\n";

my $rows = {};

while(my $row = <IN>) {
  chomp($row);
  my @f = split(/\t/,$row); 
  my $id = $f[0];
  my $chr = $f[1];
  my $start = $f[2];
  my $end = $f[3];
  my $strand = $f[4];
  my $length = $f[5];
  my $tr_name = $f[7];
  my $tr_strand = $f[11];

  my @ass = ($id,$chr,$start,$end,$strand,$length,$tr_name,$tr_strand);
  $rows->{$id.$tr_name} = \@ass;
}

foreach my $key(keys %$rows) {
  print OUT join("\t",@{$rows->{$key}})."\n";
}
