#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $fasta = $ARGV[0];

my $in = Bio::SeqIO->new(-file => $fasta,
                         -format => 'qual');

my $out = Bio::SeqIO->new(-file => ">new_$fasta",
                          -format => 'qual');

my $c = 1;

while(my $seq = $in->next_seq) {
  $seq->id($c);
  $out->write_seq($seq);
  $c++;
}
