#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $not_m2_id = 'NOT_IN_FM2.txt';
my $all_tr_fa = 'Fracy1_GeneModels_FilteredModels1_nt.fasta';
my $not_m2_fa = 'FCtr_not_m2.fa';
my $m1_m2_fa = 'FCtr_m1_m2.fa';

my $idhref = {};
open(NOT,$not_m2_id);
while(my $id = <NOT>) {
  chomp($id);
  $idhref->{$id} ++;
}

my $not_m2 = Bio::SeqIO->new(-file => ">$not_m2_fa",
                             -format => 'fasta');

my $m1_m2 = Bio::SeqIO->new(-file => ">$m1_m2_fa",
                            -format => 'fasta');

my $all_tr = Bio::SeqIO->new(-file => $all_tr_fa,
                             -format => 'fasta');
while(my $seq = $all_tr->next_seq) {
  my $id = $seq->id;
  my @id = split(/\|/,$id);
  my $pop = pop @id;
  $seq->id($pop);
  if (exists $idhref->{$seq->id}) {
    $not_m2->write_seq($seq);
    $idhref->{$seq->id} = 'OK';
  }
  else {
    $m1_m2->write_seq($seq);
  }
}

foreach my $key (keys %$idhref) {
  print "MISSING: \-$key\-\n" unless $idhref->{$key} eq 'OK';
}
