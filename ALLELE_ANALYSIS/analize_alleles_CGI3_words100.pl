#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use CGI qw/:standard/;
use Bio::SeqIO;
use Bio::AlignIO;

my $allio = 'FC_alleles.txt';
my $seqio = 'FC_promoter_1000_500.fa';
#mkdir('FC_alleles_analysis_CGI3');
chdir('FC_alleles_analysis_CGI3');
#system("cp ../$allio .");
#system("cp ../$seqio .");

my $href = build_href($allio);

my $fasta = Bio::SeqIO->new(-file => $seqio,
                            -format => 'fasta');

attach_seq($href,$fasta);

###open(OUT,">alleles.html");
open(TAB,">alleles.tab");
###print OUT start_html();

###open(TSS100,">tss100.tab");
###open(TSS500,">tss500.tab");
open(WORD100,">word100.tab");

my $WORD = {};

foreach my $key(sort{$a<=>$b}(keys%{$href->{ass}})) {

  my $id1 = $href->{ass}->{$key}->{id1};
  my $id2 = $href->{ass}->{$key}->{id2};
  next if ($href->{list}->{$id1} =~ /^\d+$/ || $href->{list}->{$id2} =~ /^\d+$/);

  my $file = "Group_$key\.fa";
###  print TSS100 "Group_$key\t";
###  print TSS500 "Group_$key\t";

  my $seq1 = Bio::Seq->new(-id => $id1,
                           -seq => $href->{list}->{$id1});

  my $seq2 = Bio::Seq->new(-id => $id2,
                           -seq => $href->{list}->{$id2});

  my $alnfile = "$file\.aln";
  $alnfile =~ s/.fa//;

###  my $prom = get_tr_aln($seq1,$seq2,"pre_$file",900,1000);
#  scan_aln("pre_$alnfile",10,10,$file,'pre',500);
  scan_aln_words("pre_$alnfile",3,1,$file,'pre',100);
#  tss100("pre_$alnfile",'pre');
#  tss500("pre_$alnfile",'pre');
###  my $tran = get_tr_aln($seq1,$seq2,"post_$file",1000,1100);
#  scan_aln("post_$alnfile",10,10,$file,'post',500);
  scan_aln_words("post_$alnfile",3,1,$file,'post',100);
#  tss100("post_$alnfile",'post');
#  tss500("post_$alnfile",'post');

###  $file =~ s/.fa//;

###  print OUT h1($file),
###    h2('Pre TSS Alignment -500/1'),
###    pre($prom),
###    h2('Post TSS Alignment 1/+500'),
###    pre($tran),
###    hr,
###    br;
}

###print OUT end_html;
###close(OUT);

print Dumper $WORD;

foreach my $pos(keys %$WORD) {
  foreach my $word(keys %{$WORD->{$pos}}) {
    print WORD100 "$pos\t$word\t".$WORD->{$pos}->{$word}->{conserved}."\t".$WORD->{$pos}->{$word}->{notconserved}."\n";
###          (($WORD->{$pos}->{$word}->{notconserved} / $WORD->{$pos}->{$word}->{conserved}) * 100)."\n";
  }
}

sub tss100 {
  my $alnf = shift;
  my $pos = shift;
  my $in = Bio::AlignIO->new(-file => $alnf);
  my $aln = $in->next_aln;
  my $seq1 = $aln->get_seq_by_pos(1);
  my $seqname = $seq1->id;
  my $slice;
  my $prestart = $aln->column_from_residue_number($seqname,401);
  my $preend = $aln->column_from_residue_number($seqname,500);
  my $poststart = $aln->column_from_residue_number($seqname,1);
  my $postend = $aln->column_from_residue_number($seqname,100);
  if($pos eq 'pre') {$slice = $aln->slice($prestart,$preend);}
  elsif($pos eq 'post') {$slice = $aln->slice($poststart,$postend);}
  print TSS100 $slice->percentage_identity;
  if($pos eq 'pre') {print TSS100 "\t";}
  elsif($pos eq 'post') {print TSS100 "\n";}
}

sub tss500 {
  my $alnf = shift;
  my $pos = shift;
  my $in = Bio::AlignIO->new(-file => $alnf);
  my $aln = $in->next_aln;
  my $seq1 = $aln->get_seq_by_pos(1);
  my $seqname = $seq1->id;
  my $slice;
  my $prestart = $aln->column_from_residue_number($seqname,1);
  my $preend = $aln->column_from_residue_number($seqname,500);
  my $poststart = $aln->column_from_residue_number($seqname,1);
  my $postend = $aln->column_from_residue_number($seqname,500);
  if($pos eq 'pre') {$slice = $aln->slice($prestart,$preend);}
  elsif($pos eq 'post') {$slice = $aln->slice($poststart,$postend);}
  print TSS500 $slice->percentage_identity;
  if($pos eq 'pre') {print TSS500 "\t";}
  elsif($pos eq 'post') {print TSS500 "\n";}
}

sub scan_aln {

  my $alnf = shift;
  my $window = shift;
  my $slide = shift;
  my $name = shift;
  my $pos = shift;
  my $seq_length = shift;
  $name =~ s/\.fa//;
  $name .= "_$pos";

  my $in = Bio::AlignIO->new(-file => $alnf);
  my $aln = $in->next_aln;
  my $seq1 = $aln->get_seq_by_pos(1); print $seq1->seq;
  my $seqname = $seq1->id;

  my $start = 1;
  my $end = $window;
  my $check;
  
  while($end<=$seq_length) {

    $end = $seq_length and $check = 1 if $end >= $seq_length;

    my $sstart = $aln->column_from_residue_number($seqname,$start);
    my $send = $aln->column_from_residue_number($seqname,$end);
    my $slice = $aln->slice($sstart,$send);
    my $id = $slice->percentage_identity;
    $id =~ s/^(\d+)\.\d+$/$1/;
    my $nseq = $slice->no_sequences;
    $id = 0 if $nseq < 2;
    print TAB "$name\t$id\n";
    $start += $slide;
    $end += $slide;
  }
}

sub scan_aln_words {

  my $alnf = shift;
  my $window = shift;
  my $slide = shift;
  my $name = shift;
  my $pos = shift;
  my $seq_length = shift;
  $name =~ s/\.fa//;
  $name .= "_$pos";
 
  my $in = Bio::AlignIO->new(-file => $alnf);
  my $aln = $in->next_aln;
  my $seq1 = $aln->get_seq_by_pos(1);
  my $seqname = $seq1->id;
  my $string = $seq1->seq;
  $string =~ s/\-//g;

  my $refseq = Bio::Seq->new(-id => $seqname,
                             -seq => $string);

  die "\n\n ERROR LENGTH: ".$refseq->length." is not equal to $seq_length\n!!!" if $refseq->length > $seq_length+1;

  my $start = 1;
  my $end = $window;
  my $check;
  
  while($end<=$seq_length) {

    $end = $seq_length and $check = 1 if $end >= $seq_length;

    my $seqword = $refseq->trunc($start,$end);
    my $word = $seqword->seq;

    my $sstart = $aln->column_from_residue_number($seqname,$start);
    my $send = $aln->column_from_residue_number($seqname,$end);
    my $slice = $aln->slice($sstart,$send);
    my $nseq = $slice->no_sequences;

    if($nseq == 2) {
      my $s1 = $slice->get_seq_by_pos(1);
      my $s2 = $slice->get_seq_by_pos(2);
      if($s1->seq eq $word) {
        if($s1->seq eq $s2->seq) {
          $WORD->{$pos}->{$word}->{conserved} ++;
        }
        else {
          $WORD->{$pos}->{$word}->{notconserved} ++;
        }
      }
      else {
        $WORD->{$pos}->{$word}->{notconserved} ++;
      }
    }
    else {
      $WORD->{$pos}->{$word}->{notconserved} ++;
    }
    $start += $slide;
    $end += $slide;
  }
}

sub get_tr_aln {
  my $seq1 = shift;
  my $seq2 = shift;
  my $file = shift;
  my $tstart = shift;
  my $tend = shift;
  my $fo = Bio::SeqIO->new(-file => ">$file",
                           -format => 'fasta');
  $fo->write_seq($seq1->trunc($tstart,$tend));
  $fo->write_seq($seq2->trunc($tstart,$tend));
  $fo->close;
  system("clustalw2 $file");
  $file =~ s/.fa//;
  open(IN,"$file\.aln");
  my @file = <IN>;
  my $string = join('',@file);
  $string =~ s/CLUSTAL 2.1 multiple sequence alignment\n//;
  return $string;
}

sub build_href {
  my $file = shift;
  open(IN,$file);
  my $href = {};
  my $c = 1;
  while(my $row = <IN>) {
    next if $row =~ /^\#/;
    chomp($row);
    my @f = split("\t",$row);
    $href->{list}->{$f[0]} ++;
    $href->{list}->{$f[1]} ++;
    $href->{ass}->{$c}->{id1} = $f[0];
    $href->{ass}->{$c}->{id2} = $f[1];
    $c ++;
  }
  return $href;
}

sub build_href_old {
  my $file = shift;
  open(IN,$file);
  my $href = {};
  while(my $row = <IN>) {
    next if $row =~ /^\#/;
    chomp($row);
    my @f = split("\t",$row);
    $href->{list}->{$f[1]} ++;
    $href->{list}->{$f[2]} ++;
    $href->{ass}->{$f[0]}->{id1} = $f[1];
    $href->{ass}->{$f[0]}->{id2} = $f[2];
  }
  return $href;
}

sub attach_seq {
  my $href = shift;
  my $fasta = shift;
  while(my $seq = $fasta->next_seq) {
    my $desc = $seq->description;
    my @f = split(/\s+/,$desc);
    my $id = $f[0];
    if(exists $href->{list}->{$id}) {
      next unless $seq->length >= 1500;
      $href->{list}->{$id} = $seq->seq;
    }
  }
}

sub attach_seq_old {
  my $href = shift;
  my $fasta = shift;
  while(my $seq = $fasta->next_seq) {
    my $desc = $seq->description;
    my @f = split(/\s+/,$desc);
    my $id = $f[0];
    if(exists $href->{list}->{$id}) {
      $href->{list}->{$id} = $seq->seq;
    }
  }
}


__END__

#!/usr/bin/perl
use strict;
use warnings;

use lib "$ENV{HOME}/src/BioPerl-1.6.1";

use Bio::AlignIO;

my $window = $ARGV[1];
my $slide = $ARGV[2];

die "ERROR:\n\tThe slide is greater than the window!\n" if $slide > $window;

my $in = Bio::AlignIO->new(-file => $ARGV[0]);

my $aln = $in->next_aln;
print $aln->length,"bp length of alignment\n";
my $start = 1;
my $end = $window;
my $check;
	
while(!$check) {

	$end = $aln->length and $check = 1 if $end > $aln->length;
	my $slice = $aln->slice($start,$end);
	my $id = $slice->percentage_identity;
	$id =~ s/^(\d+)\.\d+$/$1/;
	print "From $start to $end\t".$id.
	"% among ".$slice->no_sequences." sequences\n";
		
	$start += $slide;
	$end += $slide;
}


http://genome.lbl.gov/vista/mvista/instructions.shtml#anno

http://genome.lbl.gov/vista/mvista/submit.shtml
