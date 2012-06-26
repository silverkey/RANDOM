#!/usr/bin/perl;
use strict;
use warnings;
use Data::Dumper;

my $filename = $ARGV[0];
my $debug = $ARGV[1];
die "\n\tUSAGE: perl $0 [exonerate output] [debug]\n\n" unless $ARGV[0];
die "\n\tERROR: Cannot find the file $ARGV[0]\n\n" unless -e $ARGV[0];

open(IN,$filename);

my $row;
my $block;
my $ingff;
my $res;
my $seen;

while($row = <IN>) {
  if($row =~ /START OF GFF DUMP/) {
    $block = 1;
    $ingff = 1;
    while($block) {
      $row = <IN>;
      if($row =~ /END OF GFF DUMP/) {
        undef $ingff;
      }
      next if $row =~ /^\#/;
      if($row =~ /C4 Alignment/ || $row =~ /completed exonerate analysis/) {
        convert($res);
        print Dumper $res if $debug;
        undef $block;
        undef $ingff;
        undef $res;
      }
      $res->{gff2} .= $row if $ingff && $block;
      $res->{custom} .= $row if !$ingff && $block;
    }
  }
}

sub convert {
  my $res = shift;
  my $gff = $res->{gff2};
  my $id;
  foreach my $row(split(/\n/,$gff)) {
    my $string;
    my @field = split(/\t/,$row);
    if($field[2] eq 'gene') {
      $field[8] =~ /sequence (.+) ; gene_orientation/;
      my $preid = $1;
      $seen->{$preid} ++;
      $id = $preid.'.'.$seen->{$preid};
      $string = "$field[0]\texonerate\tmRNA\t$field[3]\t$field[4]\t$field[5]\t$field[6]\t$field[7]\tID=$id\n";
    }
    elsif($field[2] eq 'exon') {
      $string = "$field[0]\texonerate\texon\t$field[3]\t$field[4]\t$field[5]\t$field[6]\t$field[7]\tParent=$id\n";
    }
    print $string if $string;
  }
}
