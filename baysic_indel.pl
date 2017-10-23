#!/usr/bin/perl -w
#count_venn.vcf

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'fasta|f=s','capture|b=s','help|h','prefix|p=s','execdir|e=s');

unless($opt{prefix}) {
    $opt{prefix} = 'baysic';
}
my $total = 0;
if ($opt{fasta}) {
    open FASTA, "<$opt{fasta}\.fai" or die "unable to find fai index for $opt{fasta}, please index with samtools";
    while (my $line = <FASTA>) {
	chomp($line);
	my ($chr,$length,$offset,$linelen,$linebytes) = split(/\t/,$line);
	$total += $length;
    }
}else {
    open BED, "<$opt{capture}" or die "Please provide reference genome (-f) or capture bed(-b)";
    while (my $line = <BED>) {
	chomp($line);
	my ($chr,$start,$end) = split(/\t/,$line);
	$total += $end-$start+1;
    }
}

my %ct = ();
my %snpbins = ();
open CTS, ">baysic.cts" or die $!;
print CTS "Estimated sum: ".$total,"\n";

#bedtools multiinter -header -i gatk1.indel.bed platypus1.indel.bed sam1.indel.bed ssvar1.indel.bed > indel_intersect.txt

open VC, "<indel_intersect.txt" or die $!;
my $header = <VC>;
chomp($header);
my @header = split(/\t/,$header);
@callers = @header[5..$#header];
my $numcallers = scalar(@callers);
while (my $line = <VC>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my $key = join("",@row[5..$#row]);
    $ct{$key} ++;
    push @{$chrpos{$key}}, join("\t",@row[0..2]);
    $total --;
}

my @g = bits($numcallers);
$ct{$g[0]} = $total;
foreach (@g) {
    $ct{$_} = 0 unless ($ct{$_});
    print CTS join("\t",$_,$ct{$_}),"\n";
}

system("Rscript $opt{execdir}\/scripts/lca.R -c baysic.cts -s baysic.stats");

my @key1 = split(//,$g[-1]);
my @key2;
foreach $i (0..$#key1) {
    push @key2, $i if ($key1[$i] > 0);
}
open STATS, "baysic.stats" or die $!;
open BEDU, ">indels.bed" or die $!;
my @keeppos;
while (my $line = <STATS>) {
    chomp($line);
    next unless ($line =~ m/postprobs\[(\d+)\]\s+(\S+)/);
    my ($index,$prob) = ($1,$2);
    next unless ($prob >= 0.8);
    my $key = $g[$index-1];
    if ($chrpos{$key}) {
	print BEDU join("\n",@{$chrpos{$key}}),"\n";
    }
}
close BEDU;
system("sort -k 1,1 -k 2,2n indels.bed > indels.sort.bed");

sub bits { glob join "", map "{0,1}", 1..$_[0] }
