#!/usr/bin/perl -w
#parse_cnvkit_table.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'input|s=s','help|h','datadir|r=s',
			  'prefix|p=s');

open OM, "<$opt{datadir}/cancer.genelist.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}

my $file = $opt{input};

open OUT, ">$opt{prefix}.exonskip.bed" or die $!;

open IN, "<$file" or die $!;
my $header = <IN>;
chomp($header);
my @colnames = split(/\t/,$header);

while (my $line = <IN>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#row) {
	$hash{$colnames[$i]} = $row[$i];
    }
    next if ($hash{known_junction} || $hash{exons_skipped} < 1);
    next unless ($keep{$hash{genes}});
    next if ($hash{score} < 9);
    print OUT join("\t",$hash{chrom},$hash{start},$hash{end},$hash{score}),"\n";
}
close IN;
close OUT;
system(qq{bedtools intersect -wao -a $opt{prefix}.exonskip.bed -b $opt{datadir}/gencode.exon.bed > $opt{prefix}.exonskip.genes.txt});

my %skip;
open IN, "<$opt{prefix}.exonskip.genes.txt" or die $!;
while (my $line = <IN>) {
    chomp($line);
    my ($chrom,$start,$end,$readct,$chr,$a,$b,$gene,$len) = split(/\t/,$line);
    next if ($chr eq '.');
    my ($sym,$trxid,$exonnum,$strand) = split(/\|/,$gene);
    my $key = join("-",$chrom,$start,$end);
    $skip{$key}{readct} = $readct;
    $skip{$key}{gene} = $sym;
    push @{$skip{$key}{trxid}{$trxid}}, $exonnum
}

open OUT, ">$opt{prefix}.exonskip.answer.txt" or die $!;
print OUT join("\t","Gene","Chromosome","Start","End","Abberation Type","Readct","Transcripts"),"\n";

foreach my $loci (keys %skip) {
    my @trxs;
    foreach my $trxid (keys %{$skip{$loci}{trxid}}) {
	my @exonnums = sort {$a <=> $b} @{$skip{$loci}{trxid}{$trxid}};
	if (scalar(@exonnums) > 1) {
	    push @trxs, $trxid.':'.$exonnums[0].'-'.$exonnums[-1];
	}else {
	    push @trxs, $trxid.':'.$exonnums[0];
	}
    }
    print OUT join("\t",$skip{$loci}{gene},split(/-/,$loci),'Exon Skipping',
		   $skip{$loci}{readct},join("|",@trxs)),"\n";
}
