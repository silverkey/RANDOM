#!/usr/bin/perl
use strict;
use warnings;

my $datadir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/send/sel_coord';
chdir($datadir);
my @glob = glob('*coord.txt');
open(OUT,">RESULTS_COORD.xls");
my $h = 0;
foreach my $file(@glob) {
  open(IN,$file);
  my $compname = "$file";
  $compname =~ s/.coord.txt//;
  my $header = <IN>;
  chomp($header);
  unless($h) {
    print OUT "$header\tcompname\n";
    $h++;
  }
  while(my $row=<IN>) {
    chomp($row);
    print OUT "$row\t$compname\n";
  }
  close(IN);
}
close(OUT);
