#!/usr/bin/perl

# load_ucsc_tables_locally.pl
#   by Remo Sanges
#
# See POD documentation for this script at the end of the file
#

use strict;
use warnings;
use Getopt::Long; # for parsing command line arguments
use Pod::Usage;
use DBI;
use Cwd;
use Data::Dumper;

# Get command line options - if you are curious how this 
# works, check the Getopt::Long module documentation.

my($help,$tables_dir,$db,$usr,$pwd,$host,$verbose);

GetOptions('help' => \$help,
	         'tables_dir=s' => \$tables_dir,
	   			 'db=s' => \$db,
	   			 'usr=s' => \$usr,
	   			 'pwd=s' => \$pwd,
           'host=s' => \$host,
	         'verbose' => \$verbose);

if($help)  {
	pod2usage(-exitstatus=>0, -verbose=>2);
}
elsif (!$tables_dir || !$db || !$usr || !$pwd) {
	pod2usage(1);
}

my $DBH = connect_to_db($db,$usr,$pwd,$host);

my $CWD = getcwd();
chdir($tables_dir) or die "\nCannot go into $tables_dir:\n$!\n";
system('gunzip *.gz');
my @cfiles = glob('*.sql');
foreach my $sfile(@cfiles) {
	my $tfile = "$sfile";
	my $table = "$sfile";
	$tfile =~ s/\.sql$/\.txt/;
	$table =~ s/\.sql$//;
	my $create;
	my $read;
	open(IN,$sfile);
	while(my $row = <IN>) {
		if($read) {
			$create .= $row;
			$read = 0 if $row =~ /\;/;
		}
		if($row =~ /^CREATE TABLE/) {
			$create .= $row;
			$read = 1;
		}
	}
	my $csth = $DBH->prepare($create);
	$csth->execute;
	my $load = "LOAD DATA LOCAL INFILE \'$tfile\' INTO TABLE $table";
	my $lsth = $DBH->prepare($load);
	$lsth->execute;
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

load_ucsc_tables_locally.pl - create and load tables from UCSC into a MySQL DB

=head1 SYNOPSIS

perl load_ucsc_tables_locally.pl -d rs9_hg18 -u pippo -p pippo -t /home/remo/pippo

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

=item B<-t  or  --tables_dir>  <directory>
REQUIRED: Name of the directory containing both the .sql and the .txt (or .gz) for each table to load.

=back

=head1 DESCRIPTION

From the UCSC genome browser website you can download, for each table, a .sql file containing the create statements
for the table and a .txt file containing the rows to load into the table. This script will create and load those tables into a db.
The page for the human genome from which to download is http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/

=cut


