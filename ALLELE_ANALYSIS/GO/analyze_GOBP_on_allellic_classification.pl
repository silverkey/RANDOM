#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Data::Dumper;
use Bio::DB::SeqFeature::Store;

my $user = 'mysql_dev';
my $pwd = 'dEvEl0pEr';
my $GFFdbname = 'FC';
my $outfile = "alleles_GOBP.txt";
my $GOfile = '/Users/remo/ANALYSIS/MOCK/FC/ANALYSIS_MASKED/ANNOTATION/Fracy1_goinfo_FilteredModels1.tab';
my $allelefile = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/FC_alleles.txt';

my $GFFdb = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                            -user => $user,
                                            -pass => $pwd,
                                            -dsn => "dbi:mysql:$GFFdbname");

my $GO = collect_JGI_GO($GOfile);
my $idmap = get_name_protid_map($GFFdb);
my($al1,$al2) = get_gene_with_allele($allelefile);
my $noal = get_gene_with_no_allele($idmap,$al1,$al2);
prepare_table($GO,$idmap,$al1,$noal,$outfile);

sub get_name_protid_map {
  my $db = shift;
  my $href = {};
  my @gene = $db->get_features_by_type('gene');
  foreach my $gene(@gene) {
    my %p = $gene->attributes();
    my $alias = $p{Alias}[0];
    $alias =~ s/^t_//;
    my $name = $gene->name;
    $href->{$name} = $alias;
  }
  return $href;
}

sub collect_JGI_GO {
  my $tab = shift;
  my $href = {};
  open(GO,$tab) or die $!;
  while(my $row = <GO>) {
    chomp($row);
    my($pid,$goid,$goname,$gotype,$goacc) = split(/\t/,$row);
    next unless $gotype eq 'biological_process';
    $href->{protein}->{$pid}->{$goacc} ++;
    $href->{class}->{$goacc}->{count} ++;
    $href->{class}->{$goacc}->{name} = $goname;
  }
  return($href);
}

sub get_gene_with_allele {
  my $file = shift;
  my $res1 = {};
  my $res2 = {};
  open(IN,$file);
  while(my $row = <IN>) {
    chomp($row);
    my @f = split(/\t/,$row);
    my $id1 = $f[0];
    my $id2 = $f[1];
    $res1->{$id1} ++;
    $res2->{$id2} ++;
  }
  return($res1,$res2);
}

sub get_gene_with_no_allele {
  my $all = shift;
  my $al1 = shift;
  my $al2 = shift;
  my $noal = {};
  foreach my $name(keys %$all) {
    next if exists $al1->{$name};
    next if exists $al2->{$name};
    $noal->{$name} ++;
  }
  return $noal;
}

sub prepare_table {
  my $go = shift;
  my $idmap = shift;
  my $al = shift;
  my $noal = shift;
  my $outfile = shift;
  my $total_al = scalar(keys %$al);
  my $total_noal = scalar(keys %$noal);
  my $total_al_in_go = get_total_in_go($go,$al,$idmap);
  my $total_noal_in_go = get_total_in_go($go,$noal,$idmap);
  open(OUT,">$outfile");
  foreach my $goacc (keys %{$go->{class}}) {
    my $goname = $go->{class}->{$goacc}->{name};
    my $tot_cl = $go->{class}->{$goacc}->{count};
    my $al_cl = get_gene_in_class($go,$goacc,$al,$idmap);
    my $noal_cl = get_gene_in_class($go,$goacc,$noal,$idmap);
    #next unless ($al_cl + $noal_cl) > 11;
    print OUT join("\t","\"$goname\"",$total_al_in_go,$total_noal_in_go,$al_cl,$noal_cl);
    print OUT "\n";
  }
}

sub get_total_in_go {
  my $go = shift;
  my $names = shift;
  my $idmap = shift;
  my $c;
  foreach my $name(keys %$names) {
    my $id = $idmap->{$name};
    $c ++ if exists $go->{protein}->{$id};
  }
  return $c;
}

sub get_gene_in_class {
  my $go = shift;
  my $goacc = shift;
  my $names = shift;
  my $idmap = shift;
  my $c;
  foreach my $name(keys %$names) {
    my $id = $idmap->{$name};
    $c ++ if exists $go->{protein}->{$id}->{$goacc};
  }
  # ADD 1 TO AVOID PROBLEMS!
  $c++;
  return $c;
}

__END__

t = read.table(file='alleles_GO.txt',sep='\t')
colnames(t) = c('Class','All_in_GO','No_all_in_GO','All_in_class','No_all_in_class')
t = t[t[,4]>50,]
a = c()
for(i in 1:nrow(t)) {
  r = t[i,]
  pval = prop.test(as.numeric(c(r[4],r[5])),as.numeric(c(r[2],r[3])),alternative='g')$p.value
  a = c(a,pval)
}
a = p.adjust(a,'fdr')
t$adj_p = a
t[a<=0.05,]
write.table(t[a<=0.05,],file='alleles_GO_significant.xls',sep='\t',row.names=F,quote=F)
