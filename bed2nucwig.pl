# Convert a bed file to a wig file for paired end sequencing : fill in the reads
#	Extend midpoint of the reads to 30 bp on each side - assuming the central region of a nucleosome
# For calling nucleosome positions - use only nucleosome sized fragments.
die "Usage: perl bed2wig.pl <BED FILE> <WIG FILE NAME> <MIN> <MAX>\n" if(!$ARGV[3]);

$GN = 139712364; # Genome size - for making the coverage values look decent, change it to your genome size

$min=$ARGV[2];
$max=$ARGV[3];

open(LIST,$ARGV[0]) || die "INPUT $!\n";
while(chomp($l=<LIST>)){
	open(FILE,$l) || die "Bed file: $l $!\n";
	while(chomp($line=<FILE>)){
		$lno++;
		@temp = split /[\ \s\n\t]+/, $line;
		$frag_length = $temp[3]-$temp[2];
		#Assuming 6 column bed file - change if different
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$line\n";
		}elsif($frag_length>=$min && $frag_length<=$max){
			$temp[0]=~s/chr//;
			$mp = int( ($temp[1] + $temp[2])/20 + 0.5 );
			$lower=$mp-3;
			$upper=$mp+3;
			for($i=$lower;$i<=$upper;$i++){
				$largest{$temp[0]} = $i if($largest{$temp[0]}<$i);
				$read{$temp[0]}{$i}++;
				$count++;
			}
		}
		print STDERR "Count:$count\n" if($count%1000000==0);
		print STDERR "Line No:$lno\n" if($lno%10000==0);
	}
}
print STDERR "Finished reading bed file\nCount=$count\n";
close(FILE);
open(OUT,">$ARGV[1]") || die "OUT $!\n";
print OUT "track type=wiggle_0\n";
foreach $i (keys (%read) ){
	print OUT "variableStep  chrom=chr$i span=10\n";
	%thash = %{$read{$i}};
	for($j=0;$j<=$largest{$i};$j++){
		if(exists $thash{$j}){
			$normval = $thash{$j}*$GN/$count ;		
			print OUT ($j*10+1)." $normval\n";
		}
	}
}
close(OUT);
