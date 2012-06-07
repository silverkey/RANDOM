#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Data::Dumper;
use Bio::DB::SeqFeature::Store;

my $user = 'mysql_dev';
my $pwd = 'dEvEl0pEr';
my $GFFdbname = 'FC';
my $outfile = "alleles_PFAM_DOMAIN.txt";
my $PFAMfile = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/PFAM/fc_pfam.txt';
my $allelefile = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/FC_alleles.txt';

my $GFFdb = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                            -user => $user,
                                            -pass => $pwd,
                                            -dsn => "dbi:mysql:$GFFdbname");

my $PFAM = collect_MOCK_PFAM($PFAMfile);
my $idmap = get_name_protid_map($GFFdb);
my($al1,$al2) = get_gene_with_allele($allelefile);
my $noal = get_gene_with_no_allele($idmap,$al1,$al2);
prepare_table($PFAM,$idmap,$al1,$noal,$outfile);

print "\n";
print "Total  :".scalar(keys %$idmap)."\n";
print "Total 1:".scalar(keys %$al1)."\n";
print "Total 2:".scalar(keys %$al2)."\n";
print "Total N:".scalar(keys %$noal)."\n";
print "\n";
print "Done\n";

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

sub collect_MOCK_PFAM {
  my $tab = shift;
  my $href = {};
  open(PFAM,$tab) or die $!;
  while(my $row = <PFAM>) {
    chomp($row);
    my($gene_id,$pfam_id,$pfam_name,$type,$clan_id,$bit_score,$evalue) = split(/\t/,$row);
    next unless $type eq 'Domain';
    $href->{protein}->{$gene_id}->{$pfam_id} ++;
    if($href->{protein}->{$gene_id}->{$pfam_id} == 1) {
      $href->{class}->{$pfam_id}->{count} ++;
    }
    $href->{class}->{$pfam_id}->{name} = $pfam_name;
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
    $res1->{$idmap->{$id1}} ++;
    $res2->{$idmap->{$id2}} ++;
  }
  return($res1,$res2);
}

sub get_gene_with_no_allele {
  my $all = shift;
  my $al1 = shift;
  my $al2 = shift;
  my $noal = {};
  foreach my $name(values %$all) {
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
#   next unless ($al_cl + $noal_cl) > 11;
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
    #my $id = $idmap->{$name};
    #$c ++ if exists $go->{protein}->{$id};
    $c ++ if exists $go->{protein}->{$name};
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
    #my $id = $idmap->{$name};
    #$c ++ if exists $go->{protein}->{$id}->{$goacc};
    $c ++ if exists $go->{protein}->{$name}->{$goacc};
  }
  # ADD 1 TO AVOID PROBLEMS!
  $c++;
  return $c;
}

__END__

t = read.table(file='alleles_PFAM.txt',sep='\t')
colnames(t) = c('Class','All_in_GO','No_all_in_GO','All_in_class','No_all_in_class')
t = t[t[,4]>10,]
a = c()
for(i in 1:nrow(t)) {
  r = t[i,]
  pval = prop.test(as.numeric(c(r[4],r[5])),as.numeric(c(r[2],r[3])),alternative='g')$p.value
  a = c(a,pval)
}
t$p_val = a
a = p.adjust(a,'fdr')
t$adj_p = a
t[a<=0.1,]
write.table(t[a<=0.05,],file='alleles_PFAM_significant.xls',sep='\t',row.names=F,quote=F)
