#!/usr/bin/perl 
#migrate_db.pl

my $headerfile = shift @ARGV;
my @vcffiles = @ARGV;

open HEADER, "<$headerfile" or die $!;
open OUT, ">int.vcf" or die $!;
while (my $line = <HEADER>) {
  chomp($line);
  print OUT $line,"\n";;
}
close HEADER;
my @sampleorder;

my %headerlines;
foreach $vcf (@vcffiles) {
  $caller = (split(/\./,$vcf))[-3];
  open VCF, "gunzip -c $vcf|" or die $!;
  my @sampleids;
  while (my $line = <VCF>) {
    chomp($line);
    if ($line =~ m/#/) {
      if ($line =~ m/#CHROM/) {
	($chromhd, $posd,$idhd,$refhd,$althd,$scorehd,
	 $filterhd,$annothd,$formathd,@sampleids) = split(/\t/, $line);
	foreach $j (0..$#sampleids) {
	    $sampleids[$j] = (split(/\./,$sampleids[$j]))[0];
	}
	unless (@sampleorder) {
	  @sampleorder = @sampleids;
	  print OUT join("\t",$chromhd, $posd,$idhd,$refhd,$althd,$scorehd,
			 $filterhd,$annothd,$formathd,@sampleids),"\n";
	}
	next;
      }
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
    my %newgts;
    my $missingGT = 0;
  FG:foreach my $i (0..$#gts) {
      my $allele_info = $gts[$i];
      my $sid = $sampleids[$i];
      my @gtinfo = split(/:/,$allele_info);
      my %gtdata;
      if ($allele_info eq '.') {
	$newgts{$sid} = '.:.:.:.:.';
	$missingGT ++;
	next FG;
      }
      foreach my $i (0..$#deschead) {
	$gtdata{$deschead[$i]} = $gtinfo[$i];
      }
      if ($gtdata{GT} eq './.') {
	$newgts{$sid} = '.:.:.:.:.';
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
	$gtdata{AD} = join(",",$gtdata{RO},$gtdata{AO})
      } elsif (exists $gtdata{AO} && exists $gtdata{RO}) {
	$gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
	$gtdata{DP} = $gtdata{RO};
	foreach (split(',',$gtdata{AO})) {
	  $gtdata{DP} += $_;
	}
      }
      if ($gtdata{DP} && $gtdata{DP} < 3) {
	$missingGT ++;
      }
      $newgts{$sid} =  join(":",$gtdata{GT},$gtdata{DP},$gtdata{AD},$gtdata{AO},$gtdata{RO});
    }
    next if ($missingGT == scalar(@gts));
    my @newgts;
    foreach $id (@sampleorder) {
      push @newgts, $newgts{$id};
    }
    $lines{$chrom}{$pos}{$caller} = join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,@newgts),"\n";
  }
  close VCF;
}
my @callers = ('ssvar','platypus','sam','gatk','strelka2','hotspot');
if (grep(/mutect/,@vcffiles)) {
    @callers = ('sssom','pmutect','shimmer','strelka2','varscan','virmid');
}
 F1:foreach $chr (sort {$a cmp $b} keys %lines) {
   F2:foreach $pos (sort {$a <=> $b} keys %{$lines{$chr}}) {
       my $callset = join(",",keys %{$lines{$chr}{$pos}});
     F3:foreach $caller (@callers) {
	 if ($lines{$chr}{$pos}{$caller}) {
	     my ($chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,
		 $format,@gts) = split(/\t/,$lines{$chr}{$pos}{$caller});
	     $annot = $annot.";CallSet=".$callset;
	     print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
			    $filter,$annot,$format,@gts),"\n";
	     last F3;
	 }
     }
   }
}
close OUT;
