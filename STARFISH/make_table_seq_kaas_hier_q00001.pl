#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

my $kaasko = $ARGV[0];
my $fasta = $ARGV[1];
my $usage = "\n\tperl $0 [hier file 1] [fasta]\n\n";
die $usage unless scalar(@ARGV) == 2;
die $usage unless -e $kaasko;
die $usage unless -e $fasta;

print "transcript_id\tko_id\tdefinition\tsequence\n";

my $anno;
open(IN,"$kaasko");
while(my $row = <IN>) {
  if($row =~/^D\s+(comp.+)\;\s\<a href\=\"\/dbget\-bin\/www_bget\?ko\:(.+)\"\>.+\<\/a\>\s\s(.+)$/) {
    $anno->{$1} = "\t$2\t$3\t";
  }
}

my $seqio = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

while(my $seq = $seqio->next_seq) {
  my $tid = $seq->id;
  my $string = $seq->seq;
  if(exists $anno->{$tid}) {
    print $tid;
    print $anno->{$tid};
    print "$string\n";
  }
}
