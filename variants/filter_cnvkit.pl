#!/usr/bin/perl -w
#parse_cnvkit_table.pl

my $refdir = '/project/shared/bicf_workflow_ref/GRCh38/';
open OM, "<$refdir\/panel1385.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $keep{$line} = 1;
}

open ENT_ENS, "</project/shared/bicf_workflow_ref/gene2ensembl.human.txt" or die $!;
my %entrez;
my $ent_header = <ENT_ENS>;
while (my $line = <ENT_ENS>){
  chomp $line;
  my @row = split(/\t/, $line);
  $entrez{$row[2]}=$row[1];
}
close ENT_ENS;
open ENT_SYM, "</project/shared/bicf_workflow_ref/gene_info.human.txt" or die $!;
$ent_header = <ENT_SYM>;
while (my $line = <ENT_SYM>){
  chomp $line;
  my @row = split(/\t/, $line);
  $entrez{$row[2]}=$row[1];
}
close ENT_SYM;

my $file = shift @ARGV;
my $prefix = (split(/\./,(split(/\//,$file))[0]))[0];

open OUT, ">$prefix\.cnvcalls.txt" or die $!;
open BIO, ">$prefix\.data_cna_cbioportal.txt" or die $!;
print OUT join("\t","Gene","Chromosome","Start","End","Abberation Type","CN","Score"),"\n";
print BIO join("\t","Hugo_Symbol","Entrez_Gene_Id",$prefix),"\n";

open IN, "<$file" or die $!;
my $header = <IN>;
while (my $line = <IN>) {
    chomp($line);
    my ($chr,$start,$end,$geneids,$log2,$cn,$depth,
	$probes,$weight) = split(/\t/,$line);
    my %genes;
    my @ids = split(/;|,/,$geneids);
    foreach my $gid (@ids) {
	my ($key,$value) = split(/=/,$gid);
	if ($key eq 'ensembl_gn' || $key eq 'identifier') {
	    $genes{$value} = 1 if $keep{$value};
	}
    }
    my $len = sprintf("%.1f",($end-$start)/1000);
    next if ($cn == 2) || scalar(keys %genes) < 1;
    my $abtype = 'amplification';
    $abtype = 'loss' if ($cn < 2);
    foreach $gene (keys %genes) {
	$cn_cbio = $cn -2;
	$cn_cbio = 2 if ($cn > 4);
	print BIO join("\t",$gene,$entrez{$gene},$cn_cbio),"\n";
	print OUT join("\t",$gene,$chr,$start,$end,$abtype,$cn,$weight),"\n";
    }
}
close IN;
close OUT;
close BIO;
