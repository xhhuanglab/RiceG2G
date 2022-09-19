#perl Select_MLM_Peak.pl SH_Leaf_angle
use POSIX qw(log log10);
use List::Util qw/max/;
$input="./data/farmCPU/$ARGV[0].MLM.csv ";
open IN,"$input" or die "can't open $input";
while (<IN>){
	chomp();
	if($.>1){
		@line =split/\,/,$_;
		if($line[7]<=0.0000000167  and $line[7] ne "NA"){##按阈值筛选marker
			$line[0]=~s/\"//g;
			push @snp,$line[0]; 
			$chromosome{$line[0]}=$line[1];
			$position{$line[0]}=$line[2];
			$Pvalue{$line[0]}=-log10($line[7]);
			$seq{$line[0]}=$#snp;
		}
	}
}
close IN;
$chr=0;$bin=0;
foreach $i(0..$#snp){
	if($chromosome{$snp[$i]} == $chr){
		if($position{$snp[$i]}-$pos <=1000000){
			push @{"bin$bin"},$snp;
			$pos=$position{$snp[$i]};
			$snp=$snp[$i];
		}else{
			push @{"bin$bin"},$snp;
			$pos=$position{$snp[$i]};
			$snp=$snp[$i];
			$bin++;
		}
			
			
	}else{
		push @{"bin$bin"},$snp;
		$chr=$chromosome{$snp[$i]};
		$pos=$position{$snp[$i]};
		$snp=$snp[$i];
		$bin++;
	}
}
if($chromosome{$snp[$#snp]} == $chr){
	if($position{$snp[$#snp]}-$position{$snp[$#snp-1]} <=1000000){
		push @{"bin$bin"},$snp[$#snp];
	}else{
		push @{"bin$bin"},$snp[$#snp];
	}
}else{
	push @{"bin$bin"},$snp[$#snp];
}

$num=0;
##输出每个区段最高的peak
open OUT,">$ARGV[0].MLM.peak";
foreach $i(1..$bin){
	if($#{"bin$i"}==0){
		$key=${"bin$i"}[0];
		$key=${"bin$i"}[0];
		push @peak,$key;
		print 	OUT "$key	$chromosome{$key}	$position{$key}	$Pvalue{$key}\n";
	}else{
		@P=();
		foreach $j(@{"bin$i"}){
			push @P,$Pvalue{$j};
		}
		$max_p= max @P;
		foreach $j(@{"bin$i"}){
			if($Pvalue{$j} eq $max_p){
				push @peak,$j;
				print 	OUT "$j	$chromosome{$j}	$position{$j}	$Pvalue{$j}\n";
				last;
			}
		}
	}
}
 
		
						
							