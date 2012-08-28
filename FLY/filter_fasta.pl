#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "\n\tUSAGE: perl $0 [fasta] [list]\n\n";
my $fasta = $ARGV[0];
my $list = $ARGV[1];
die $usage unless scalar(@ARGV) == 2;
die $usage unless -e $fasta;
die $usage unless -e $list;

my $outname = "$fasta";
$outname =~ s/\..+$//;
$outname .= "_$list";
$outname =~ s/\..+$/.fasta/;

my $seqin = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

my $seqout = Bio::SeqIO->new(-file => ">$outname",
                            -format => 'fasta');

open(IN,$list);

my $href = {};
while(my $row = <IN>) {
  chomp($row);
  $href->{$row} ++;
}

while(my $seq = $seqin->next_seq) {
  if(exists $href->{$seq->id}) {
    $seqout->write_seq($seq);
  }
}

