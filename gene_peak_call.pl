#!/usr/bin/perl

die "Usage: perl gene_peak_call.pl <WIG> <TSS> <OUT-FILE> <window>\n" if(!$ARGV[3]);

BEGIN { push @INC, 'perl_library' }
use ngs;

$tread = &ngs::readwig($ARGV[0]);
%read  = %{$tread};

$window=$ARGV[3];

# read TSSs
open(NUC,$ARGV[1]) || die "PEAK $!\n";
while($i=<NUC>){ 
	chomp($i);
	@temp=split/[\t\n\ ]+/,$i;
	$peak= (int($temp[3]/10+0.5))*10+1;
	#print "$temp[1]\t$temp[2]\t$peak\n";
	$chr=$temp[1];
	$chr=~s/chr//;
	$quartile=$temp[0];
	$str=$temp[4];
	$#tarray=-1;
	$#time=-1;
	for($m=$peak-$window;$m<=$peak+$window;$m+=10){
		$val=0;
		if(exists $read{$chr}{$m}){
			$val=$read{$chr}{$m};
		}
		push(@tarray,$val);
		push(@time,$m);
	}
	$star = &smooth_array(\@tarray,2);
	@starray = @{$star};
	$tx = &peaks(\@starray,\@time);
	@qarray=@{$tx};
	print "track type=wiggle_0\nvariableStep  chrom=chr$chr\n";
	for($k=0;$k<=$#qarray;$k++){
		print "$qarray[$k] 45\n";
	}

}
sub smooth_array {
	#open(RA,">ra") || die "$!\n";
	#print STDERR "smooth start\n";
	my @tar = @{$_[0]};
	my $rwin = $_[1];
	my @rar;
	for(my $ii = 0; $ii<=$#tar ; $ii++){
		if($ii>=$rwin && $ii <= $#tar-$rwin){
			my $rval=0; my $rct=0;
			for(my $jj=$ii-$rwin; $jj<=$ii+$rwin;$jj++){
				$rval+=$tar[$jj]; $rct++;
			}
			$rval /= $rct;
			$rar[$ii] = $rval;
		}else{
			$rar[$ii]=$tar[$ii];
		}
		#print RA "$time[$ii] $rar[$ii] $tar[$ii]\n";
	}
	#print STDERR "smooth end\n";
	return(\@rar);
	#close(RA);
}

sub peaks {
  my @ar = @{$_[0]};
	my @ti = @{$_[1]};
	my $sum=0;
	$#der=-1; $#p=-1; $#runAve=-1;$#runSD=-1;
	#open(DER, ">der") || die "$!\n";
	for(my $ii=0;$ii<=$#ar;$ii++){ # Populate derivative array
		if($ii<3 || $ii > $#ar-2){
			$derivative=10;
		}else{
			$h=($ti[$ii+2] - $ti[$ii-2]);
			$h*=-1 if($h<0);
			$derivative = ($ar[$ii-2] - 8*$ar[$ii-1] + 8*$ar[$ii+1] - $ar[$ii+2])/(3*$h); 		
		}
		push(@der,$derivative);
		#print DER "$ti[$ii]\t$der[$ii]\n";
	}
	close(DER);
	#print STDERR "Size of derivative Array $#der\n";
	for(my $ii=0;$ii<=$#ar;$ii++){ # Populate running ave/running SD array
		my $s=0;$no=0;
		for( my $jj=$ii-10;$jj<$ii+10;$jj++){			
			if($jj>=0 && $jj<=$#ar){
				$s+= $ar[$jj];
				$no++;
			}		
		}
		$s/=$no;
		push(@runAve,$s);
		$n=0;
		for( my $jj=$ii-10;$jj<$ii+10;$jj++){			
			$n+= ($ar[$jj]-$s)*($ar[$jj]-$s) if($jj>=0 && $jj<=$#ar); 		
		}
		$n = sqrt($n/$no);
		push(@runSD,$n);
	}
	# derivative changes from positive to negative - peak
	$left_of_peak=0;
	$#pe=-1; $#p=-1;my $cutoff;
	for($ii=50;$ii<=$#der-50;$ii++){
		$cutoff =  $runAve[$ii] + 0.5*$runSD[$ii];
		#print STDERR "$time[$ii] $cutoff\n";
		#print STDERR "$ii Derivative= $der[$ii] ; Val = $ar[$ii] ; Running Ave = $runAve[$ii] ; Running SD = $runSD[$ii] ; Cut-off = ".($runAve[$ii]+0.5*$runSD[$ii])."\n";
		if($der[$ii]>0.04){
			$left_of_peak=1;
			$#p=-1;
			push(@p,$time[$ii]);
		}elsif($left_of_peak==1 && $ar[$ii-1] > ($runAve[$ii]+0.5*$runSD[$ii]) && $der[$ii]>-0.04 && $der[$ii]<0.04){
			push(@p,$time[$ii]);
		}elsif($der[$ii]<=-0.04 && $left_of_peak==1){
			push(@p,$time[$ii]);
			$left_of_peak=0;
			if($#p>-1){
				my $s=0;
				for(my $jj=0;$jj<=$#p;$jj++){
					$s+=$p[$jj];
				}
				$s=10*int( $s/(($#p+1)*10) + 0.5);
				push(@pe,$s);
				$#p=-1;
			}
		}
		#print STDERR "Left-of-peak $left_of_peak\n";
	}
	return(\@pe);
}
