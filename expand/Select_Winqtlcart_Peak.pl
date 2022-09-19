#perl Select_Winqtlcart_Peak.pl SH_Heading_date
use List::Util qw/max/;
if($ARGV[0]=~/SH_(\S+)/){
	$trait=$1;
	$area="SH";
}elsif($ARGV[0]=~/HN_(\S+)/){
	$trait=$1;
	$area="HN";
}elsif($ARGV[0]=~/HZ_(\S+)/){
	$trait=$1;
	$area="HZ";
}
open OUT,">$ARGV[0].winqtlcart.peak";
foreach $NAM(1..15){
	@nam_lod=();@nam_lod_tr=();@lod=();@chr_marker=();$hang_star=();$id_pheno=();$hang_star=0;
	$file='./data/WinQTLcart/'."$area".'_Nam'."$NAM".'-C.qrt';
	open IN,"$file" or die "can't open $file";
	while (<IN>){
		chomp();
		$hang++;$array=$_;
		if($array=~/-trait\s+(\d+)\s+Analyzed trait \[(\S+)\]/){
			$trait_id=$1;$trait_name=$2;
			if($trait_name eq $trait){
				$hang_star++;
				$id_pheno=$trait_id;
			}
			if(defined($id_pheno) and $hang_star>1){
				last;
			}
		}elsif($array=~/^(\s+)?(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+/ and $hang_star>=1 ){
			$hang_star++;
			$lod= $5*2.5/11.5; 
			if(defined ($lod[$2][$3])){
				if($lod >$lod[$2][$3]){
					$lod[$2][$3]=$lod;
				}
			}else{
				$lod[$2][$3]=$lod;
			}
			$chr_marker[$2]=$3;
		}
		
	}
	close IN;
	$chr=0;@position=();@position_star=();@position_end=();@marker_num=();$num=0;
	open IN1,"./data/WinQTLcart/N$NAM.bin.list.map";
	while (<IN1>){
		chomp();
		$line=$_;
		if($line=~/chr(\d+)	(\S+)/){
			if($1==$chr){
				$num++;
				$position[$chr][$num]=$2;
				$position_star[$chr][$num]=$position;
				$position_end[$chr][$num]=$2*2-$position;
				$position=$position_end[$chr][$num];
			}else{
				$marker_num[$NAM][$chr]=$num;
				$chr++;$num=1;$position=0;
				$position[$chr][$num]=$2;
				$position_star[$chr][$num]=$position;
				$position_end[$chr][$num]=$2*2;
				$position=$position_end[$chr][$num];
			}
		}
	}
	close IN;
	$marker_num[$NAM][$chr]=$num;
	foreach $chr(1..12){
		$star=0;
		foreach $i(1..$marker_num[$NAM][$chr]){
			if(defined($lod[$chr][$i])){
				$star=$lod[$chr][$i];
			}else{
				$lod[$chr][$i]=$star;
			}
		}
	}
	@snp_chr=();@snp_win=();
	foreach $chr(1..12){
		foreach $i(1..$marker_num[$NAM][$chr]){
			if($lod[$chr][$i] >= 3.5){
				push @snp_chr,$chr; 
				push @snp_win,$i; 
			}
		}
	}
	###确定peak
	$chr=0;
	foreach $i(0..$bin){
		@{"bin$i"}=();
	}
	$bin=0;
	foreach $i(0..$#snp_chr){
		if( $snp_chr[$i]  == $chr){
			if($position[$snp_chr[$i]][$snp_win[$i]]-$pos <=10){
				push @{"bin$bin"},$i;
				$pos=$position[$snp_chr[$i]][$snp_win[$i]];
				$snp=$i;
			}else{
				$pos=$position[$snp_chr[$i]][$snp_win[$i]];
				$snp=$i;
				$bin++;
			}
		}else{
			push @{"bin$bin"},$snp;
			$chr=$snp_chr[$i];
			$pos=$position[$snp_chr[$i]][$snp_win[$i]];
			$snp=$i;
			$bin++;
		}
	}
	push @{"bin$bin"},$#snp;
	$num=0;
	##输出每个区段最高的peak
	foreach $i(1..$bin){
		if($#{"bin$i"}==0){
			$key=${"bin$i"}[0];
			$pos_star=$position_star[$snp_chr[$key]][$snp_win[$key]]/10;
		    $pos_end=$position_end[$snp_chr[$key]][$snp_win[$key]]/10;
			print 	OUT "NAM$NAM	$snp_chr[$key]	$pos_star	$pos_end	$lod[$snp_chr[$key]][$snp_win[$key]]\n";
		}else{
			@P=();
			foreach $j(@{"bin$i"}){
				push @P,$lod[$snp_chr[$j]][$snp_win[$j]];
			}
			$max_p= max @P;
			foreach $j(@{"bin$i"}){
				if($lod[$snp_chr[$j]][$snp_win[$j]] eq $max_p){
					$key=$j;
					 $pos_star=$position_star[$snp_chr[$key]][$snp_win[$key]]/10;
					 $pos_end=$position_end[$snp_chr[$key]][$snp_win[$key]]/10;
					 
					print OUT "NAM$NAM	$snp_chr[$key]	$pos_star	$pos_end	$lod[$snp_chr[$key]][$snp_win[$key]]\n";
					last;
				}
			}
		}
	}
}
 



