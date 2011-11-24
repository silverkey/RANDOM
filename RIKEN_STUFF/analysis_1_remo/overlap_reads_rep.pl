#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long; # for parsing command line arguments
use Pod::Usage;
use DBI;
use Data::Dumper;
use Bio::Range;

# Get command line options - if you are curious how this 
# works, check the Getopt::Long module documentation.

my($help,$coord_file,$outfile,$db,$usr,$pwd,$host,$verbose);

GetOptions('help' => \$help,
           'coord_file=s' => \$coord_file,
           'outfile=s' => \$outfile,
           'db=s' => \$db,
           'usr=s' => \$usr,
           'pwd=s' => \$pwd,
           'host=s' => \$host,
           'verbose' => \$verbose);

if($help)  {
  pod2usage(-exitstatus=>0, -verbose=>2);
}
elsif (!$coord_file || !$outfile || !$db || !$usr || !$pwd) {
  pod2usage(1);
}

my $DBH = connect_to_db($db,$usr,$pwd,$host);

open(IN,$coord_file);
my $INheader = <IN>;

open(OUT,">$outfile");
my @OUTheader = qw(id chr start end strand length 
                   rep.hit.length overlap.length rep.rna.ratio rna.rep.ratio
                   rep.chr rep.start rep.end rep.strand 
                   rep.class rep.family rep.name
                   rep.in.start rep.in.end rep.in.left);
print OUT join("\t",@OUTheader);
print OUT "\n";

while(my $element = <IN>) {

  chomp($element);
  my($id,$chr,$start0,$end,$strand) = split(/\t/,$element);
  my $overlap = get_overlap($element);

  next unless $overlap;

  foreach my $or(@$overlap) {

    my @ovrow = ($id, $chr, $start0, $end, $strand, ($end-$start0),
                 $or->{rep_length}, $or->{overlap_length}, $or->{rep_rna_ratio}, $or->{rna_rep_ratio},
                 $or->{genoName}, $or->{genoStart}, $or->{genoEnd}, $or->{strand},
                 $or->{repClass}, $or->{repFamily}, $or->{repName},
                 $or->{repStart}, $or->{repEnd}, $or->{repLeft});

    print OUT join("\t",@ovrow);
    print OUT "\n";
  }
}

sub get_overlap {

  my $element = shift;
  my($id,$chr,$start0,$end,$strand) = split(/\t/,$element);
  my $res = [];

  my $sth = $DBH->prepare("SELECT * FROM $chr\_rmsk WHERE ".
                          "genoStart <= $end AND genoEnd >= $start0");
  $sth->execute;

  while(my $row = $sth->fetchrow_hashref) {

    my $r1 = Bio::Range->new(-start => $start0,
                             -end => $end);

    my $r2 = Bio::Range->new(-start => $row->{genoStart},
                             -end => $row->{genoEnd});

    my $intersection = $r1->intersection($r2);

    # Subtract 1 to the intersection to solve the 0 based start UCSC / 1 based start BIOPERL
    $row->{overlap_length} = $intersection->length - 1;
    $row->{rep_length} = $row->{genoEnd}-$row->{genoStart};
    $row->{rep_rna_ratio} = sprintf("%.2f", $row->{overlap_length} / $row->{rep_length} * 100);
    $row->{rna_rep_ratio} = sprintf("%.2f", $row->{overlap_length} / ($end-$start0) * 100);

    # Here we can put some kind of analysis about the repeat
    # $row->{rep_info} = check_rep($row);

    push(@$res,$row);
  }
  return $res;
}

sub connect_to_db {
  my $db = shift;
  my $usr = shift;
  my $pwd = shift;
  my $host = shift;
  my $dsn = 'dbi:mysql:'.$db;
  $dsn .= ':'.$host if $host; # IN THE CURRENT DBI POD VERSION THERE IS THE '@' IN THE PLACE OF ':'
  my $dbh = DBI->connect($dsn,$usr,$pwd,{PrintWarn=>1,PrintError=>1,RaiseError=>1}) or die $DBI::errstr;
  return $dbh;
}

__END__

=head1 NAME

overlap_reads_rep.pl - makes the overlap between a coord.txt file from OSC2tab.R and the repeats tables from UCSC database

=head1 SYNOPSIS

perl overlap_reads_rep.pl -d hg18 -u pippo -c pippo -c /home/pippo/PARVALBUMIN/coord.txt -o rep.overlap.txt

=head1 OPTIONS

=over 8

=item B<-d  or  --db>  <databse name>
REQUIRED

=item B<-u  or  --usr>  <databse user>
REQUIRED

=item B<-p  or  --pwd>  <databse password>
REQUIRED

=item B<-ho  or  --host>  <databse host>
OPTIONAL
      
=item B<-c  or  --coord_file>  <filename>
REQUIRED: name of the coord file against which to overlap repeats
    
=item B<-o  or  --outfile>  <filename>
REQUIRED: name of the file to create
    
=back

=head1 DESCRIPTION
   
Some deeper description here ..................

=cut

