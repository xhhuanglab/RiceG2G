#perl SelectFarmCPUPeakByCor.pl Plant_height
use Statistics::Basic qw(:all);
use Statistics::RankCorrelation;
use Statistics::Distributions;
use Statistics::LSNoHistory;
use Statistics::Descriptive;
use POSIX qw(log log10);
use List::Util qw/max/;
$input="/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/rMVP/unmiss/$ARGV[0].FarmCPU.csv";
open IN,"$input" or die "can't open $input";
while (<IN>){
	chomp();
	if($.>1){
		@line =split/\,/,$_;
		if($line[7]<=0.0000000167 and $line[7] ne "NA"){##按阈值筛选marker 
		#if($line[7]<=0.00001 and $line[7] ne "NA"){
			$line[0]=~s/\"//g;
			push @snp,$line[0]; 
			$chromosome{$line[0]}=$line[1];
			$position{$line[0]}=$line[2];
			$Pvalue{$line[0]}=sprintf('%.02f',-log10($line[7]));
			if($Pvalue{$line[0]} eq "inf"){
				$Pvalue{$line[0]}=309.00;
			}
			$seq{$line[0]}=$#snp;
			$effect{$line[0]}=$line[5];
			$se{$line[0]}=$line[6];
		}
	}
}
close IN;
$chr=0;$bin=0;
foreach $i(0..$#snp){
	if($chromosome{$snp[$i]} == $chr){
		if($position{$snp[$i]}-$pos <=2000000){
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
	if($position{$snp[$#snp]}-$position{$snp[$#snp-1]} <=2000000){
		push @{"bin$bin"},$snp[$#snp];
	}else{
		push @{"bin$bin"},$snp[$#snp];
	}
}else{
	push @{"bin$bin"},$snp[$#snp];
}



open IN3,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/rMVP/unmiss/Parent_unmiss.tped" or die "can't open Parent_unmiss.tped";
while (<IN3>){
	chomp();
	if($_=~/\d+	(snp\d+)	0	\d+/){
		$key=$1;
		if(defined($Pvalue{$key})){
			@array=split/\s+/,$_;
			foreach $i(0..(($#array+1)/2)-3){
				if($array[$i*2+4] eq $array[$i*2+5]){
					if($array[$i*2+4] == 0){
						$parent_geno{$key}.=0;
					}elsif($array[$i*2+4] == 1){
						$parent_geno{$key}.=1;
					}elsif($array[$i*2+4] == 2){
						$parent_geno{$key}.=2;
					}else{
						$parent_geno{$key}.="error1";
					}
				}else{
					$parent_geno{$key}.="error2";
				}
			}
		}
	}
}
close IN3;

##将每个区段内按window大小划分，输出最高的peak
foreach $i(1..$bin){
	if($#{"bin$i"}==0){
		push @{"select_bin$i"},${"bin$i"}[0];
	}else{
		@current_P=();@current_snp=();
		$num=1;
		$key0=${"bin$i"}[0];
		foreach $j(0..$#{"bin$i"}){
			$key=${"bin$i"}[$j];
			if($position{$key}-$position{$key0}<=2000000*$num){
				push @current_P,$Pvalue{$key};
				push @current_snp,$key;
			}else{
				$maxp=max @current_P;
				foreach $m(0..$#current_P){
					if($current_P[$m] eq $maxp){
						push @{"select_bin$i"},$current_snp[$m];
						last;
					}
				}
				@current_P=();@current_snp=();
				$num=int(($position{$key}-$position{$key0})/2000000)+1;
				push @current_P,$Pvalue{$key};
				push @current_snp,$key;
			}
		}
		$maxp=max @current_P;
		foreach $m(0..$#current_P){
			if($current_P[$m] eq $maxp){
				push @{"select_bin$i"},$current_snp[$m];
				last;
			}
		}
	}
}
foreach $i(1..$bin){
	if($#{"select_bin$i"}==0){
		$key=${"select_bin$i"}[0];
		$marker[$marker_num][0]=$i;
		$marker[$marker_num][1]=$key;
		$marker[$marker_num][2]=$chromosome{$key};
		$marker[$marker_num][3]=$position{$key};
		$marker[$marker_num][4]=$Pvalue{$key};
		$marker[$marker_num][5]=$parent_geno{$key};
		$marker_num++;
	}else{
		foreach $j(@{"select_bin$i"}){
			$marker[$marker_num][0]=$i;
			$marker[$marker_num][1]=$j;
			$marker[$marker_num][2]=$chromosome{$j};
			$marker[$marker_num][3]=$position{$j};
			$marker[$marker_num][4]=$Pvalue{$j};
			$marker[$marker_num][5]=$parent_geno{$j};
			$marker_num++;
		}
	}
}

@sort_marker=sort{$a->[0]<=>$b->[0] or $a->[2]<=>$b->[2] or $b->[4]<=>$a->[4] or $a->[3]<=>$b->[3]}@marker;

$star=1;
foreach $i(0..$marker_num-1){
	if($sort_marker[$i][0] == $star){
		push @{"sort_bin$star"},$sort_marker[$i][1];
	}else{
		$star++;
		push @{"sort_bin$star"},$sort_marker[$i][1];
	}
}
$final_num=0;
foreach $i(1..$star){
	$binnum=$#{"sort_bin$i"};
	if($binnum==0){
		$snpName=${"sort_bin$i"}[$binnum];
		$final[$final_num][0]=$chromosome{$snpName};
		$final[$final_num][1]=$position{$snpName};
		$final[$final_num][2]=$Pvalue{$snpName};
		$final[$final_num][3]=$parent_geno{$snpName};
		$final[$final_num][4]=$effect{$snpName};
		$final[$final_num][5]=$se{$snpName};
		push @select_key,"$final[$final_num][0]	$final[$final_num][1]";
		$final_num++; 
	}else{
		@name1=();
		$m=1;
		foreach $j(0..$binnum){
			${"name$m"}[$j]=${"sort_bin$i"}[$j];
		}
		
		foreach $j(1..$binnum+1){
			if( $#{"name$j"} > -1){
				if($#{"name$j"}>0){
					$name1=${"name$j"}[0];
					$snpName=${"sort_bin$i"}[$binnum];
					$final[$final_num][0]=$chromosome{$name1};
					$final[$final_num][1]=$position{$name1};
					$final[$final_num][2]=$Pvalue{$name1};
					$final[$final_num][3]=$parent_geno{$name1};
					$final[$final_num][4]=$effect{$name1};
					$final[$final_num][5]=$se{$name1};
					push @select_key,"$final[$final_num][0]	$final[$final_num][1]";
					$final_num++; 
					$a=$j+1;
					@{"name$a"}=();
					foreach $m(1..$#{"name$j"}){
						$name2=${"name$j"}[$m];
						$correlation=&cor($parent_geno{$name1},$parent_geno{$name2});
						if($correlation<=0.2){
							push @{"name$a"},$name2;
						}
					}
				}else{
					$name1=${"name$j"}[0];
					$final[$final_num][0]=$chromosome{$name1};
					$final[$final_num][1]=$position{$name1};
					$final[$final_num][2]=$Pvalue{$name1};
					$final[$final_num][3]=$parent_geno{$name1};
					$final[$final_num][4]=$effect{$name1};
					$final[$final_num][5]=$se{$name1};
					push @select_key,"$final[$final_num][0]	$final[$final_num][1]";
					$final_num++; 
					@{"name$j"}=();
				}
			}else{
				last;
			}
		}
	}
}
@sort_final=sort{$a->[0]<=>$b->[0] or $a->[1]<=>$b->[1]}@final;
foreach $i(1..12){
	open IN4,"/data/home/mjchen/rice/NAM_denovo_assemble/mummer4/indel_blast/chr$i-nam-Effect.snp";
	while (<IN4>){
		chomp();
		if($_=~/Chr(\d+)	(\d+)	(\S+)	(\S+)	/){
			$key="$1	$2";
			$ref{$key}=$3;
			$alt{$key}=$4;
		}
	}
	close IN4;
}

open IN5,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/rMVP/unmiss/slelct_LD_NAM_Magic_0.8.tped" or die "can't open the file slelct_LD_NAM_Magic_0.8.tped";
while (<IN5>){
	chomp();
	
	if($_=~/(\d+)	snp\d+	\d+	(\d+)	/){
		$snpid="$1	$2";
		$ref=0;$alt=0;$wrong=0;
		if(grep {$snpid eq $_} @select_key){
			@tped=split/\s+/,$_;
			foreach $i(4..$#tped){
				if($tped[$i] == 1){
					$ref++;
				}elsif($tped[$i] == 2){
					$alt++;
				}else{
					$wrong++;
				}
			}
			$ref_frequency{$snpid}=$ref/2;
			$alt_frequency{$snpid}=$alt/2;
		}
	}
}
close IN5;

	

open OUT,">$ARGV[0].FarmCPU.peak" or die "can't open the output file";
print OUT "Chromosome	Position	-log10(Pvalue)	Parent_geno	1_genotype	2_genotype	1_frequency	2_frequency	Effect	SE\n";
foreach (0..$final_num-1){
	$key="$sort_final[$_][0]	$sort_final[$_][1]";
	print OUT "$sort_final[$_][0]	$sort_final[$_][1]	$sort_final[$_][2]	$sort_final[$_][3]	$ref{$key}	$alt{$key}	$ref_frequency{$key}	$alt_frequency{$key}	$sort_final[$_][4]	$sort_final[$_][5]\n"
}
close OUT;


sub cor {
	@parent1=split//,$_[0];
	@parent2=split//,$_[1];
	$cor=Statistics::RankCorrelation-> new( \@parent1, \@parent2 );
	my $stat_snp=Statistics::Descriptive::Full->new();
	my $stat_corr=Statistics::Descriptive::Full->new();
	$stat_snp->add_data(@parent1);
	$stat_corr->add_data(@parent2);
	$standard_snp=$stat_snp-> standard_deviation;
	$standard_corr=$stat_corr-> standard_deviation;
	$df=$cor->{'size'} - 2;
	$reg = Statistics::LSNoHistory->new;#计算pearson相关系数
	$reg->append_arrays([@parent1],[@parent2]);#计算pearson相关系数
	$rho=$reg->pearson_r;#计算pearson相关系数
	$r2=$rho*$rho;
	
}



