#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;

# DEFAULTS
my $DIR = '/Users/remo/ANALYSIS/mmmir124';
my $UTR_FASTA = '/Users/remo/ANALYSIS/mmmir124/mm_pans.txt';
my $MIR_FASTA = '/Users/remo/ANALYSIS/mmmir124/mus_musculus.mirbase';
my $OUT = 'scan_mmpans_mm_mirbase.xls';

my $result = GetOptions ('utr_fasta|u=s' => \$UTR_FASTA,
                         'mir_fasta|m=s' => \$MIR_FASTA,
                         'output|o=s' => \$OUT,
                         'dir|d=s' => \$DIR);

chdir($DIR);
open(OUT,">$OUT");
print OUT "gene_id\ttranscript_id\tmiRNA_id\t7mer\-A1\t7mer-m8\t8mer\tTOT\n";

my $UTR = Bio::SeqIO->new(-file => $UTR_FASTA,
                          -format => 'fasta');

my $MIR = Bio::SeqIO->new(-file => $MIR_FASTA,
                          -format => 'fasta');

my $mirna_seqs = get_mirna();

while(my $seq = $UTR->next_seq) {
  my @ref = split(/\|/,$seq->id);
  my $gid = $ref[0]; # gene stable id
  my $tid = $ref[1]; # transcript stable id
  foreach my $id(keys %$mirna_seqs) {
    my $counts = scan_utr($seq->seq,$mirna_seqs->{$id});
    my $sum = $counts->{'7mer-A1'}->{count}+$counts->{'7mer-m8'}->{count}+$counts->{'8mer'}->{count};
    print OUT join("\t",($gid,$tid,$id,$counts->{'7mer-A1'}->{count},$counts->{'7mer-m8'}->{count},$counts->{'8mer'}->{count},$sum));
    print OUT "\n";
  }
}

close(OUT);

exit;

sub get_mirna {
  my $href = {};
  while(my $seq = $MIR->next_seq) {
    $href->{$seq->id} = $seq->seq;
  }
  return $href;
}

sub scan_utr {

	my $utr = shift;
	my $mirna = shift;
	my $motif = '';
	my $motifs = {};

	# CHANGE 'U' WITH 'T'
	$mirna =~ tr/uU/tT/;

	# EXTRACT NUCLEOTIDES 2..8 FROM MIRNA
	my @bases = split(//,$mirna);
	my $core = join('',@bases[1..7]);
	my $seed = join('',@bases[1..6]);

	my $core_rev = reverse($core);
	$core_rev =~ tr/ACGTacgt/TGCAtgca/;
	my $seed_rev = reverse($seed);
	$seed_rev =~ tr/ACGTacgt/TGCAtgca/;
	
	# 7mer-A1
	# NUCLEOTIDES 2..7 INVERTED + 'A' AT THE END
	$motifs->{'7mer-A1'}->{count} = 0;
	$motif = "$seed_rev".'A';
	$motifs->{'7mer-A1'}->{motif} = $motif;
	while($utr =~ /$motif/ig) {
		$motifs->{'7mer-A1'}->{count} ++;
	}

	# 7mer-m8
	# NUCLEOTIDES 2..8 INVERTED
	$motifs->{'7mer-m8'}->{count} = 0;
	$motif = "$core_rev";
	$motifs->{'7mer-m8'}->{motif} = $motif;
	while($utr =~ /$motif/ig) {
		$motifs->{'7mer-m8'}->{count} ++;
	}

	# 8mer
	# NUCLEOTIDES 2..8 INVERTED + 'A' AT THE END
	$motifs->{'8mer'}->{count} = 0;
	$motif = "$core_rev".'A';
	$motifs->{'8mer'}->{motif} = $motif;
	while($utr =~ /$motif/ig) {
		$motifs->{'8mer'}->{count} ++;
	}

	for(my $c=1; $c<=$motifs->{'8mer'}->{count}; $c++) {
		$motifs->{'7mer-A1'}->{count} --;
		$motifs->{'7mer-m8'}->{count} --;
	}

	return $motifs;
}
