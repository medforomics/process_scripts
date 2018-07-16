#!/usr/bin/perl 
#migrate_db.pl

my $vcf = shift @ARGV;
my $out = $vcf;
$out =~ s/\.vcf/.gt.vcf/g;

open VCF, "<$vcf" or die $!;
open OUT, ">$out" or die $!;
while (my $line = <VCF>) {
    chomp($line);
    $line =~ s/ID:/ID=/g;
    if ($line =~ m/#CHROM/) {
	print OUT join("\t",$line,'FORMAT','NORMAL','TUMOR'),"\n";
    }elsif ($line =~ m/#/) {
	print OUT $line,"\n";
    }else {
	my ($chrom, $pos,$id,$ref,$alt,$score,
	    $filter,$annot) = split(/\t/, $line);
	foreach $a (split(/;/,$annot)) {
	    my ($key,$val) = split(/=/,$a);
	    $hash{$key} = $val;
	}
	$normalgt=join(":",'0/0',$hash{NDP},$hash{NAC},$hash{NDP}-$hash{NAC},
		       join(',',$hash{NDP}-$hash{NAC},$hash{NAC}));
	$tumorgt=join(":",'0/0',$hash{DDP},$hash{DAC},$hash{DDP}-$hash{DAC},
		       join(',',$hash{DDP}-$hash{DAC},$hash{DAC}));
	print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,
		       'GT:DP:AO:RO:AD',$normalgt,$tumorgt),"\n";
    }
}
