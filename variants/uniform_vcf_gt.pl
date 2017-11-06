#!/usr/bin/perl 
#migrate_db.pl

my $vcf = shift @ARGV;
my $outfile = $vcf;
$outfile =~ s/vcf.gz/uniform.vcf/;
open OUT, ">$outfile" or die $!;
open VCF, "gunzip -c $vcf|" or die $!;
while (my $line = <VCF>) {
    chomp($line);
    if ($line =~ m/#/) {
	print OUT $line,"\n";
	next;
    }
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
	my ($key,$val) = split(/=/,$a);
	$hash{$key} = $val;
	}
    my @deschead = split(/:/,$format);
    my $newformat = 'GT:DP:AD:AO:RO';
    my @newgts = ();
    my $missingGT = 0;
  FG:foreach my $allele_info (@gts) {
      my @gtinfo = split(/:/,$allele_info);
      my %gtdata;
      if ($allele_info eq '.') {
	  push @newgts, '.:.:.:.:.';
	  $missingGT ++;
	  next FG;
      }
      foreach my $i (0..$#deschead) {
	  $gtdata{$deschead[$i]} = $gtinfo[$i];
      }
      if ($gtdata{DP} == 0 || $gtdata{GT} eq './.') {
	  push @newgts, '.:.:.:.:.';
	  $missingGT ++;
	  next FG;
      }
      if ($gtdata{AD}){
	  ($gtdata{RO},@alts) = split(/,/,$gtdata{AD});
	  $gtdata{AO} = join(",",@alts);
	  $gtdata{DP} = $gtdata{RO};
	  foreach (@alts) {
	      $gtdata{DP} += $_;
	  }
      } elsif (exists $gtdata{NR} && exists $gtdata{NV}) {
	      $gtdata{DP} = $gtdata{NR}; 	
	      $gtdata{AO} = $gtdata{NV};
	      $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
      } elsif (exists $gtdata{AO} && exists $gtdata{RO}) {
	  $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
	  $gtdata{DP} = $gtdata{RO};
	  foreach (split(',',$gtdata{AO})) {
	      $gtdata{DP} += $_;
	  }
      }
      if ($gtdata{DP} && $gtdata{DP} < 5) {
	  $missingGT ++;
      }
      push @newgts, join(":",$gtdata{GT},$gtdata{DP},$gtdata{AD},$gtdata{AO},$gtdata{RO});
  }
    next if ($missingGT == scalar(@gts));
    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,@newgts),"\n";
}
close VCF;
