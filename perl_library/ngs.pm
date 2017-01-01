package ngs;
use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = 1.00;

sub writeWig{
	my ($href,$outfile)=@_;
	my %Mseq  = %{$href};
	open(OUT,">$outfile") || die "OUT writeWig $outfile $!\n";
	print OUT "track type=wiggle_0\n";
	foreach my $i (keys (%Mseq) ){
		print OUT "variableStep  chrom=chr$i span=10\n";
		my %thash = %{$Mseq{$i}};
		foreach my $j ( sort {$a<=>$b} keys(%thash) ){
			print OUT "$j $Mseq{$i}{$j}\n";
		}
	}
	close(OUT);
	return(1);
}

sub readYwig{
#wig format: Chr Coord Value
	my @ytemp; my %Ywig; my $ywig_line;	
	open(YWIG,$_[0]) || die "YWIG $!\n";
	while($ywig_line=<YWIG>){
		@ytemp = split /[\ \s\n\t]+/, $ywig_line;
		$ytemp[0]=~s/chr//;
		$Ywig{$ytemp[0]}{$ytemp[1]}=$ytemp[2]; # $Ywig{Chromosome}{coordinate} = value
	}
	close(YWIG);
	return (\%Ywig);
}

sub subwig{
# outwig = subwig(peakfile, window, wig_file)
	my @temp; my %wig; my $wig_line; my $chrid; my %pos;
	my $window = $_[1];
	print "$window\n";
	open(FILE,$_[0]) || die "Peak file $!\n";
	while(chomp(my $peak_line=<FILE>)){
		@temp = split /[\ \s\n\t]+/, $peak_line;
		$chrid = $temp[1];
		$chrid=~s/chr//;
		my $peak = (int($temp[3]/10+0.5))*10+1;
		print "$chrid $peak\n";
		for(my $i=$peak-$window;$i<=$peak+$window;$i+=10){
			$pos{$chrid}{$i} = 1;
		}
	}
	close(FILE);

	open(WIG,$_[2]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//;
			$chrid=$temp[2];
			print "$chrid\n";
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1] if(exists $pos{$chrid}{$temp[0]}); # $wig{Chromosome}{coordinate} = value
		}
	}
	close(WIG);
	return (\%wig);
}

sub make_peakhash {
	# hash_ref = make_peakhash(peak_file, window)
	my $chrid; my %pos; my $window = $_[1]; my $peak_file = $_[0];
	my @temp;
	open(FILE,$peak_file) || die "Peak file $!\n";
	while(chomp(my $peak_line=<FILE>)){
		@temp = split /[\ \s\n\t]+/, $peak_line;
		$chrid = $temp[1];
		$chrid=~s/chr//;
		my $peak = int(($temp[3]/10) + 0.5)*10 + 1;
		print "$chrid $peak\n";
		for(my $i=$peak-$window;$i<=$peak+$window;$i+=10){
			$pos{$chrid}{$i} = 1;
		}
	}
	print STDERR "Peak Filler Filled\n";
	close(FILE);
	return(\%pos);

}

sub subpairs{
# jnk = subpairs(peak_hash, pairs_file, out_file)
	my $href=$_[0];
	my %pos  = %{$href};
	my $pair_file = $_[1];
	my $outfile=$_[2];
	my @temp;
	print STDERR "pairs: $pair_file out_pairs: $outfile\n";
	open(OUT,'>>',$outfile) || die "OUT writeWig $outfile $!\n";
	open(PAIRS,$pair_file) || die "WIG $!\n";
	while(my $bed_line=<PAIRS>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		my $st = int(($temp[2]/10) + 0.5)*10 + 1;
		my $en = int(($temp[3]/10) + 0.5)*10 + 1;
		$temp[0]=~s/chr//;
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif(exists $pos{$temp[0]}{$st} || exists $pos{$temp[0]}{$en}){
			print OUT $bed_line;
		}
	}
	close(PAIRS);
	close(OUT);
	print STDERR "Done with pairs file\n";
	return (1);
}



sub readDwig{
#wig format: Chr Coord Value
	my @ytemp; my %Dwig; my $dwig_line;	
	open(DWIG,$_[0]) || die "YWIG $!\n";
	while($dwig_line=<DWIG>){
		@ytemp = split /[\ \s\n\t]+/, $dwig_line;
		$ytemp[0]=~s/chr//;
		$Dwig{$ytemp[0]}{$ytemp[1]}=$ytemp[3] if($ytemp[3]!=0); # $Dwig{Chromosome}{coordinate} = value
	}
	close(DWIG);
	return (\%Dwig);
}

sub readbedgraph{
#format: chr start_coord end_coord value
	my @temp; my %bedgraph; my $bed_line; my $chrid;
	open(BEDGRAPH,$_[0]) || die "WIG $!\n";
	while(my $bed_line=<BEDGRAPH>){
		@temp = split /[=\ \s\n\t]+/, $bed_line;
		if($bed_line=~/track/){
		}else{
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			my $st=( int( ($temp[1]+1) / 10 ) )*10+1;
			my $en=( int(  $temp[2]    / 10 ) )*10+1;
			for(my $i=$st;$i<=$en;$i+=10){
				$bedgraph{$chrid}{$i}=$temp[3];
			}
		}
	}
	close(BEDGRAPH);
	return (\%bedgraph);
}

sub readpairs_fillten{
	my ($pair_file,$min,$max)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0;
	my $i; my $j;
	open(PAIRS,$pair_file) || die "WIG $!\n";
	while(my $bed_line=<PAIRS>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif($temp[5]>=$min && $temp[5]<=$max){
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			my $lower=int(($temp[2]/10) + 0.5)*10 + 1; 
			my $upper=int(($temp[3]/10) + 0.5)*10 + 1;
			for($i=$lower;$i<=$upper;$i+=10){
				$twig{$chrid}{$i}++;
				$ct++;
			}
		}
	}
	print STDERR "Lines used: $ct\n";
	my $GN = 139712364; # Genome size
	#normalize
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			my $normval = $twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	close(PAIRS);
	return (\%wig);
}


sub readwig{
#wig format: Coord Value ; chr in the header only
	my @temp; my %wig; my $wig_line; my $chrid;
	open(WIG,$_[0]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//g;
			$chrid=$temp[2];
			print STDERR "$chrid\n";
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1]; # $Ywig{Chromosome}{coordinate} = value
		}
	}
	close(WIG);
	return (\%wig);
}

sub readwigArray{
#wig format: Coord Value ; chr in the header only
	my @temp; my %wig; my $wig_line; my $chrid; my @wig_ar;
	open(WIG,$_[0]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//;
			$chrid=$temp[2];
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1]; # $Ywig{Chromosome}{coordinate} = value
			push(@wig_ar,$temp[1]);
		}
	}
	close(WIG);
	return (\%wig,\@wig_ar);
}

sub readFlank{
#flank format:
# 0 - chr
# 1 - left flank  (genomic, not based on direction of transcription)
# 2 - right flank ( " ")
# 3 - FBgn name
# 4 - filler
# 5 - strand
# 6 - start (TSS)
# 7 - end (TES)
# 8 - flank size (optional)
	my @ftemp; my %gene; my $fline;
	open(FLANKFILE,$_[0]);
	while(chomp($fline=<FLANKFILE>)){
		@ftemp = split /[\ \s\n\t]+/, $fline;
		$gene{$ftemp[0]}{$ftemp[3]}{'start'}  = $ftemp[6];
		$gene{$ftemp[0]}{$ftemp[3]}{'end'}    = $ftemp[7];
		$gene{$ftemp[0]}{$ftemp[3]}{'strand'} = $ftemp[5];
		#convert genomic direction into transcription direction for flanks
		if($ftemp[5] eq "+"){
			$gene{$ftemp[0]}{$ftemp[3]}{'lf'} = $ftemp[1];
			$gene{$ftemp[0]}{$ftemp[3]}{'rf'} = $ftemp[2];
		}
		if($ftemp[5] eq "-"){
			$gene{$ftemp[0]}{$ftemp[3]}{'lf'} = $ftemp[2];
			$gene{$ftemp[0]}{$ftemp[3]}{'rf'} = $ftemp[1];
		}
	}
	close(FLANKFILE);
	return(\%gene);
}

sub normTSS_bp{
# Genehash, wighash, left flank size, right flank size
	my ($href1,$href2,$Lflank,$Rflank)=@_;
	my %thash = %{$href1};
	my %Mseq  = %{$href2};
	my $lst; my $rst;
	foreach my $i ( keys(%thash) ) {
		foreach my $j ( keys(%{$thash{$i}}) ){
			if($thash{$i}{$j}{'strand'} eq "+"){
				$lst=$thash{$i}{$j}{'start'} + 25;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}+$Rflank if($rst - $thash{$i}{$j}{'start'} > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m<=$rst;$m++){
					$count++;
					my $val = 0;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} +\n";
			}elsif($thash{$i}{$j}{'strand'} eq "-"){
				$lst=$thash{$i}{$j}{'start'}-25;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}-$Rflank if($thash{$i}{$j}{'start'}-$rst > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m>=$rst;$m--){
					my $val = 0;
					$count++;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} -\n";
			}
		}
	}
	return(\%thash);
}


sub normGB_bp{
# Genehash, wighash, Start_Pos, Flank size
	my ($href1,$href2,$Lst,$Rflank)=@_;
	my %thash = %{$href1};
	my %Mseq  = %{$href2};
	my $st; my $rst;
	foreach my $i ( keys(%thash) ) {
		foreach my $j ( keys(%{$thash{$i}}) ){
			if($thash{$i}{$j}{'strand'} eq "+"){
				$st=$thash{$i}{$j}{'start'} + $Lst;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}+$Rflank if($rst - $thash{$i}{$j}{'start'} > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$st;$m<=$rst;$m++){
					my $val = 0;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
						$count++;
						$sum+=$val;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} +\n";
			}elsif($thash{$i}{$j}{'strand'} eq "-"){
				$st=$thash{$i}{$j}{'start'}-$Lst;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}-$Rflank if($thash{$i}{$j}{'start'}-$rst > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$st;$m>=$rst;$m--){
					my $val = 0;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
						$count++;
						$sum+=$val;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} -\n";
			}
		}
	}
	return(\%thash);
}


sub normTSS{
# Genehash, wighash, left flank size, right flank size
	my ($href1,$href2,$Lflank,$Rflank)=@_;
	my %thash = %{$href1};
	my %Mseq  = %{$href2};
	my $lst; my $rst;
	foreach my $i ( keys(%thash) ) {
		foreach my $j ( keys(%{$thash{$i}}) ){
			if($thash{$i}{$j}{'strand'} eq "+"){
				$lst=(int($thash{$i}{$j}{'lf'}/10))*10+1;
				$lst=(int( ($thash{$i}{$j}{'start'}-$Lflank)/10))*10+1 if($thash{$i}{$j}{'start'}-$lst > $Lflank);
				$rst=(int($thash{$i}{$j}{'rf'}/10))*10+1;
				$rst=(int( ($thash{$i}{$j}{'start'}+$Rflank)/10))*10+1 if($rst - $thash{$i}{$j}{'start'} > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m<=$rst;$m+=10){
					if(exists $Mseq{$i}{$m}){
						$sum+=$Mseq{$i}{$m};
						$count++;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'}\n";
			}elsif($thash{$i}{$j}{'strand'} eq "-"){
				$lst=(int($thash{$i}{$j}{'lf'}/10))*10+1;
				$lst=(int( ($thash{$i}{$j}{'start'}+$Lflank)/10))*10+1 if($lst - $thash{$i}{$j}{'start'} > $Lflank);
				$rst=(int($thash{$i}{$j}{'rf'}/10))*10+1;
				$rst=(int( ($thash{$i}{$j}{'start'}-2000)/10))*10+1 if($thash{$i}{$j}{'start'}-$rst > 2000);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m>=$rst;$m-=10){
					if(exists $Mseq{$i}{$m}){
						$sum+=$Mseq{$i}{$m};
						$count++;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'}\n";
			}
		}
	}
	return(\%thash);
}

sub median {
  my @ar     = sort { $a <=> $b } @{$_[0]};
  my $mp     = int ( $#ar / 2);
  my $median = $ar[$mp];
	my $num = $#ar+1;
  return ($median,$num);
}



sub mean_se {
  my @ar = @{$_[0]};
  my $sum=0;
	my $x;
  foreach my $ii (@ar){
    $sum+=$ii;
  }
  my $mean=$sum/($#ar+1);
  $sum=0;
  foreach my $ii (@ar){
    $x=($ii-$mean)**2;
    $sum+=$x;
  }
  my $se=(sqrt($sum))/($#ar+1);
	my $num = $#ar+1;
  return ($mean,$se,$num);
}

sub mean_sd {
  my @ar = @{$_[0]};
  my $sum=0;
	my $x;
  foreach my $ii (@ar){
    $sum+=$ii;
  }
  my $mean=$sum/($#ar+1);
  $sum=0;
  foreach my $ii (@ar){
    $x=($ii-$mean)**2;
    $sum+=$x;
  }
  my $se=sqrt($sum/($#ar+1));
	my $num = $#ar+1;
  return ($mean,$se,$num);
}


1;
