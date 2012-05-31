#!/usr/bin/perl
use strict;
use warnings;

# SET UP THE VARIABLES CONTAINING THE NAMES OF THE FOLDERS AND THE LOCATIONS OF THE FILES

# directory containing the fasta and the quality files - HERE WILL BE SAVED THE RESULTS ALSO
my $fasta_dir = '/media/LOCAL_DATA/RESEARCH/ANALYSIS/OCTOPUS/COMBINED_ASSEMBLY/ASSEMBLY/CAP3';

# name of the fasta file - THE QUALITY IS TO BE IN THE SAME DIR WITH EXACTLY THE SAME NAME OF THE FASTA PLUS THE .QUAL EXTENSION
my $fasta = 'octopus_transcripts_out.unpadded.fasta';

# go into directory with the fasta
chdir($fasta_dir);

# launch the analysis
system("cap3 $fasta \> cap.out");
