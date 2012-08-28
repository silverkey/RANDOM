#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# The script will prepare a comprehensive gff3 file to be
# loaded into a gbrowse. It expect to find in the folder
# in which it runs:
#
# - fasta files containing the sequences of the genome
# - gff3 files containing the results of exonerate/repeatmasker etc.
# - tab delimited text file containing annotations such as the
#   output of annocript
#
# After the execution you simply have to load into MySQL the fasta
# file with the genomic regions and the generated gff3
#
# Example:
# perl bp_seqfeature_load.PLS -c -f -a DBI::mysql -d [DATABASE NAME] [GENOME FASTA] [GENOME GFF3]
#
# To create the database in MySQL you have to run:
# - create database [DATABASE NAME];
# - grant all privileges on volvox.* to [YOUR USERNAME]@localhost;
# - grant select on [DATABASE NAME].* to nobody@localhost;
#
# After having installed gbrose it is suggested to change some privileges:
# - chown [YOUR USERNAME] /var/lib/gbrowse/databases/databases
# - chown [YOUR USERNAME] /etc/gbrowse
#
# Then to create the configuration file for a new database:
# - cd /var/lib/gbrowse/databases/databases
# - mkdir [DATABASE NAME]
# - chmod go+rwx [DATABASE NAME]
#
# In the folder created you need to put the genomic fasta and the gff3
# created with this script. Then you need to edit the /etc/gbrowse/GBrowse.conf
# file and create a /etc/gbrowse/[DATABASE NAME].conf file.

# ----------------------------------- #
# BEGINNING OF HARD CODED PARAMETERS! #
# ----------------------------------- #

# Name of the gff3 output to create
my $OUTPUT = 'pm_v1.0_10000_all.gff3';

# Index of the element containing the transcript id
# in the annotation file
my $IDFIELD = 1;

# Index of the element containing the annotation
# to collect in the annotation file
my $ANNOFIELD = 4;

# String corresponding to the lack of annotations
# in the annotation file
my $NULLANNO = 'NULL';

# ----------------------------- #
# END OF HARD CODED PARAMETERS! #
# ----------------------------- #

system("rm $OUTPUT");
die "\nCANNOT REMOVE THE OUTPUT\n" if -e $OUTPUT;

# All the fasta files into the folder
my @fa = glob('*.fa');
my @fasta = glob('*.fasta');

# All the gff3 files into the folder
my @gff3 = glob('*.gff3');

# Possible annotation data from tables
# We assume they are all tab delimited text
my @anno = glob('*.xls');
my @tsv = glob('*.tsv');
my @txt = glob('*.txt');

push(@fasta,@fa);
push(@anno,@tsv);

print STDERR "\n";
print STDERR scalar(@fasta)." fasta       : @fasta\n";
print STDERR scalar(@gff3). " gff3        : @gff3\n";
print STDERR scalar(@anno). " annotations : @anno\n\n";

# Currently we can only collect a single annotation for each transcript
my $href = {};
foreach my $anno (@anno) {
  open(IN,$anno);
  # Get rid of the header (we assume the annotation file contains it)
  my $head = <IN>;
  while(my $row = <IN>) {
    my @field = split(/\t/,$row);
    next if $field[$ANNOFIELD] eq $NULLANNO;
    $href->{$field[$IDFIELD]} = $field[$ANNOFIELD];
  }
}
  
open(OUT,">$OUTPUT");

# The fasta should contains the sequence of the genome.
# Here we create the 'contig' lines of the gff3 corresponding
# to the entire contigs/scaffolds/chromosomes.
foreach my $fasta (@fasta) {
  my $seqio = Bio::SeqIO->new(-file => $fasta,
                              -format => 'fasta');
  while(my $seq = $seqio->next_seq) {
    my $last_field = 'ID='.$seq->id.';Name='.$seq->id;
    print OUT join("\t",$seq->id,'contig','contig',1,$seq->length,'.','.','.',$last_field);
    print OUT "\n";
  }
}

# The gff3 should contains bioinformatics annotations such as
# exonerate results, repeats masker results and so on.
# Their format must be in gff3 in order to be then able to use
# the bioperl/gbrowse scripts.
foreach my $gff3 (@gff3) {
  open(IN,$gff3);
  while(my $row = <IN>) {
    next if $row =~ /^\#/;
    chomp($row);
    my @field = split(/\t/,$row);
    # This is a bit specific to the gff3 we prepare from the
    # exonerate output. Probably when will adopt 'maker2' this
    # script will be different or maybe not needed at all.
    if($field[2] eq 'mRNA') {
      my $id = $field[8];
      $id =~ s/ID\=//;
      my $name = "$id";
      # In order to permit multiple mapping of the same EST
      # we add a .1/.2/.3.... to the id of the mapped ESTs. 
      # This extension is taken out in the 'Name' field. It also
      # not usually present into the annotation file.
      $name =~ s/\.\d+$//;
      my $note;
      if(exists $href->{$name}) {
        $note = $href->{$name};
        # Change all the '=' and ';' character from the annotation
        # because they a specific meaning in gff3 and are badly
        # interpreted by the bioperl/gbrowse scripts.
        $note =~ s/\;/\,/g;
        $note =~ s/\=/\:/g;
        print OUT "$row;Name=$name;Note=$note\n";
      }
      else {
        print OUT "$row;Name=$name\n";
      }
    }
    else {
      print OUT "$row\n";
    }
  }
}
