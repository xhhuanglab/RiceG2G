#perl WinQTLcart_lod_geno_rap.pl SH_Heading_date
use POSIX qw(log log10);
use List::Util qw/max/;
use Statistics::Basic qw(:all);
use Statistics::RankCorrelation;
use Statistics::Distributions;
use Statistics::LSNoHistory;
use Statistics::Descriptive;
$chr_length{1}=43270923;
$chr_length{2}=35937250;
$chr_length{3}=36413819;
$chr_length{4}=35502694;
$chr_length{5}=29958434;
$chr_length{6}=31248787;
$chr_length{7}=29697621;
$chr_length{8}=28443022;
$chr_length{9}=23012720;
$chr_length{10}=23207287;
$chr_length{11}=29021106;
$chr_length{12}=27531856;
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
foreach $NAM(1..15){
	@nam_lod=();@nam_lod_tr=();@lod=();@chr_marker=();$hang_star=();$id_pheno=();$hang_star=0;
	$file="./data/WinQTLcart/$area".'_Nam'."$NAM".'-C.qrt';
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
	$chr=0;
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
			foreach $j($position_star[$chr][$i]..$position_end[$chr][$i]){
				$lod_win[$NAM][$chr][$j]=$lod[$chr][$i];
				
			}
		}
	}
}

			

$chr=1;$star=0;
open IN2,"./basic_data/Rice_IRGSP-1.0.gff3" or die "can't open Rice_IRGSP-1.0.gff3";
while (<IN2>){
	chomp();
	if($_=~/chr(\d+)	irgsp1_locus	gene	(\d+)	(\d+)	/){
		$chromosome=$1+1-1;$pos=int(($2+$3)/2);
		$geno_star=$2;
		$geno_end=$3;
		if($_=~/ID=(\S+);Name=(\S+);/){
			$geno=$1;
			push @geno,$geno;
			push @geno_star,$geno_star;
			push @geno_end,$geno_end;
			push @geno_chr,$chromosome;
			if($chromosome==$chr){
				if(int($pos/100000)==$star){
					foreach $NAM(1..15){
						$key="$geno-$NAM";
						$lod{$key}=$lod_win[$NAM][$chr][$star];
					}
				}else{
					$star=int($pos/100000);
					foreach $NAM(1..15){
						$key="$geno-$NAM";
						$lod{$key}=$lod_win[$NAM][$chr][$star];
					}
				}
			}else{
				$chr++;$star=0;
				if(int($pos/100000)==$star){
					foreach $NAM(1..15){
						$key="$geno-$NAM";
						$lod{$key}=$lod_win[$NAM][$chr][$star];
					}
				}else{
					$star=int($pos/100000);
					foreach $NAM(1..15){
						$key="$geno-$NAM";
						$lod{$key}=$lod_win[$NAM][$chr][$star];
					}
				}
			}
				
		}
	}
}
close IN2;
open OUT,">$ARGV[0]_NAM_LOD_geno_rap.info";

print OUT "Geno_id	Chr	Geno_star	Geno_end	NAM1	NAM2	NAM3	NAM4	NAM5	NAM6	NAM7	NAM8	NAM9	NAM10	NAM11	NAM12	NAM13	NAM14	NAM15\n";
foreach $i(0..$#geno){
	print OUT "$geno[$i]	$geno_chr[$i]	$geno_star[$i]	$geno_end[$i]";
	foreach $j(1..15){
		$key= "$geno[$i]-$j";
		if(defined ($lod{$key})){
			print OUT "	$lod{$key}";
		}else{
			print OUT "	NA";
		}
	}
	print OUT "\n";
}





