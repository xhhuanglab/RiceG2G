#perl select_FarmCPUpeak_info_RAP.pl SH_Plant_height shoot_seeding_stage
use List::Util qw/min max/;
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
###转座子信息
open IN1,"all.locus_brief_info.7.0" or die "can't open all.locus_brief_info.7.0";
while (<IN1>){
	chomp();
	@line =split/\t+/,$_;
	if($line[0]=~/Chr(\d+)/){
		if($line[6] eq "Y"){
			$annotation{$line[1]}="$1".'transposon';
		}
	}
}
close IN1;
###看家基因信息
open IN2,"house-keeping_gene.list" or die "house-keeping_gene.list";
while (<IN2>){
	chomp();
	if($.>1){
		$house_keeping{$_}="House_keeping";
	}
}
close IN2;
###重要基因信息
open IN3,"repressed_geno.list" or die "repressed_geno.list";
while (<IN3>){
	chomp();
	if($.>1){
		$Important{$_}="Important";
	}
}
close IN3;
###拟南芥同源信息
open IN4,"/data/home/mjchen/rice/NAM_denovo_assemble/Candidate_gene_capture/Homology_alignment/Basic_Information_Arabidopsis_to_rice.txt" or die "can't open Basic_Information_Arabidopsis_to_rice.txt";
while (<IN4>){
	chomp();
	if($.>1){
		@Arabidopsis=split/\t+/,$_;
		if($_=~/LOC(\S+)	AT\S+/){
			if($Arabidopsis[1] eq ""){
				$Arabidopsis[1]="*";
			}
			if($Arabidopsis[2] eq ""){
				$Arabidopsis[2]="*";
			}
			if($Arabidopsis[3] eq ""){
				$Arabidopsis[3]="*";
			}
			if($Arabidopsis[5] eq ""){
				$Arabidopsis[5]="*";
			}
			if($Arabidopsis[6] eq ""){
				$Arabidopsis[6]="*";
			}
			$AT{$Arabidopsis[0]}="$Arabidopsis[1]=>$Arabidopsis[2]=>$Arabidopsis[3]=>$Arabidopsis[5]=>$Arabidopsis[6]";
		}
	}
}
close IN4;

open RAP4,"/data/home/mjchen/rice/NAM_denovo_assemble/Candidate_gene_capture/Homology_alignment/Basic_Information_Arabidopsis_to_rice_rap.txt" or die "can't open Basic_Information_Arabidopsis_to_rice.txt";
while (<RAP4>){
	chomp();
	if($.>1){
		@Arabidopsis=split/\t+/,$_;
		if($_=~/^Os(\S+)	AT\S+/){
			if($Arabidopsis[1] eq ""){
				$Arabidopsis[1]="*";
			}
			if($Arabidopsis[2] eq ""){
				$Arabidopsis[2]="*";
			}
			if($Arabidopsis[3] eq ""){
				$Arabidopsis[3]="*";
			}
			if($Arabidopsis[5] eq ""){
				$Arabidopsis[5]="*";
			}
			if($Arabidopsis[6] eq ""){
				$Arabidopsis[6]="*";
			}
			$AT{$Arabidopsis[0]}="$Arabidopsis[1]=>$Arabidopsis[2]=>$Arabidopsis[3]=>$Arabidopsis[5]=>$Arabidopsis[6]";
		}
	}
}
close RAP4;


###基因位置信息
open IN5,"/data/home/mjchen/rice/reference/MSU/Rice_MSUv7.gff3" or die "can't open Rice_MSUv7.gff3";
while (<IN5>){
	chomp();
	if($_=~/chr(\d+)	MSU_osa1r7	gene	(\d+)	(\d+)	/){
		$chr=$1+1-1;$star=$2;$end=$3;
		if($_=~/ID=(\S+);Name/){
			$geno=$1;
			foreach $i($star..$end){
				$position[$chr][$i]=$geno;
			}
			$geno_star{$geno}=$star;
			$geno_end{$geno}=$end;
			$geno_chr{$geno}=$chr;
		}
	}
}
close IN5;
open IN1,"/data/home/mjchen/rice/reference/MSU/RAP-MSU_2021-11-11.txt" or die "can't open pair file";
while (<IN1>){
	chomp();
	@pair=split/\t+/,$_;
	if($pair[1] eq "None"){
		push @RAP,$pair[0];
	}
}
close IN1;
open INRAP,"/data/home/mjchen/rice/reference/IRGSP-1.0_representative/locus.gff" or die "can't open Rice_MSUv7.gff3";
while (<INRAP>){
	chomp();
	if($_=~/chr(\d+)	irgsp1_locus	gene	(\d+)	(\d+)	/){
		$chr=$1+1-1;$star=$2;$end=$3;
		if($_=~/ID=(\S+);Name=(\S+);/){
			$geno=$1;
			if(grep {$geno eq $_} @RAP){
				foreach $i($star..$end){
					if(defined($position[$chr][$i])){
					}else{
						$position[$chr][$i]=$geno;
					}
				}
				$geno_star{$geno}=$star;
				$geno_end{$geno}=$end;
				$geno_chr{$geno}=$chr;
			}
		}
	}
}
close INRAP;


#加入RNA表达量数据
#open RNA,"/data/home/mjchen/rice/NAM_denovo_assemble/RNA_seq_database/Final_TPM_matrix_MSU.txt" or die "can't open RNA data";
open RNA,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/RNA_db/combineRAPvsMSU_RNAseq.txt" or die "can't open RNA data";
while (<RNA>){
	chomp();
	if($.>1){
		@rna=split/\s+/,$_;
		foreach $i(4..$#rna){
			$rna_geno{$rna[0]} .="	$rna[$i]";
		}
	}
}
close RNA;
##15个群体winqtlcart定位LOD值
open IN10,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/NAM_Winqtlcart/$ARGV[0]_NAM_LOD_geno.info" or die "can't open $ARGV[0]_NAM_LOD_geno.info";
while (<IN10>){
	chomp();
	if($.>1){
		@win=split/\s+/,$_; 
		$pos{$win[0]}="$win[1]:$win[2]-$win[3]";
		foreach $i(4..$#win){
			push @{"LOD$win[0]"},$win[$i];
		}
	}
}
close IN10;
open RAP10,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/NAM_Winqtlcart/$ARGV[0]_NAM_LOD_geno_rap.info" or die "can't open $ARGV[0]_NAM_LOD_geno.info";
while (<RAP10>){
	chomp();
	if($.>1){
		@win=split/\s+/,$_; 
		$pos{$win[0]}="$win[1]:$win[2]-$win[3]";
		foreach $i(4..$#win){
			push @{"LOD$win[0]"},$win[$i];
		}
	}
}
close RAP10;

###MLM混合模型信息
open IN11,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/rMVP/unmiss/$ARGV[0].MLM.peak";
while (<IN11>){
	chomp();
	@mlm=split /\s+/,$_;
	$num++;
	$mlm_peak[$num][0]=$mlm[1];
	$mlm_peak[$num][1]=$mlm[2];
	$mlm_peak[$num][2]=$mlm[3];
}
close IN11;
###winqtlcarpeak信息
$win_peaknum=0;
open IN12,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/NAM_Winqtlcart/$ARGV[0].winqtlcart.peak" or die "can't open INput 12";
while (<IN12>){
	chomp();
	@win_peak=split/\s+/,$_;
	$win_peaknum++;
	$win_peakNAM[$win_peaknum]=$win_peak[0];
	$win_peakchr[$win_peaknum]=$win_peak[1];
	$win_peakstar[$win_peaknum]=$win_peak[2]*1000000;
	$win_peakend[$win_peaknum]=$win_peak[3]*1000000;
	$win_peaklod[$win_peaknum]=$win_peak[4];
}
close IN12;
####RAP的注释信息
open RAPID,"/data/home/mjchen/rice/reference/MSU/RAP-MSU_2021-11-11.txt";
while (<RAPID>){
	chomp();
	@rap=split/\t+/,$_;
	@msuid=();
	@msuid=split/\,/,$rap[1];
	foreach $i(@msuid){
		if($i=~/(LOC_\S+).\d+/){
			$msu_rapid{$rap[0]}=$1;
		}
	}
}
close RAPID;
open ARPann,"/data/home/mjchen/rice/reference/IRGSP-1.0_representative/locus.gff";
while (<ARPann>){
	if($_=~/ID=(\S+);Name=/){
		$id=$1;
	}
	if($_=~/Note=(.*)/){
		$ann=$1;
	}
	$msu_rap_ann{$msu_rapid{$id}}=$ann;
	$msu_rap_ann{$id}=$ann;
}
close ARPann;
open IN9,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/rMVP/unmiss/$ARGV[0].FarmCPU.peak" or die "can't open farmCPU.info";
while (<IN9>){
	chomp();
	if($.>1){
		@farm=split/\s+/,$_;
		$mlm_info=();
		foreach $i(1..$num){
			if($mlm_peak[$i][0]==$farm[0]){
				if(abs ($mlm_peak[$i][1]-$farm[1]) <=300000){
					$mlm_info.="Chr$mlm_peak[$i][0]-$mlm_peak[$i][1]-$mlm_peak[$i][2]";
				}
			}
		}
		if($mlm_info eq ""){
			$mlm_info="*";
		}
		$mlm_info{$QTL}=$mlm_info;
		$QTL_chr{$QTL}=$farm[0];
		$win_peak_info=();
		@win_peak_star=();
		@win_peak_end=();
		@win_peak_lod=();
		@new_win_peak_star=();
		@new_win_peak_end=();
		@NAM_id=();
		foreach $i(1..$win_peaknum){
			if($farm[0]==$win_peakchr[$i]){
				if(abs($win_peakstar[$i]-$farm[1])<=500000 or abs($win_peakend[$i]-$farm[1])<=500000){
					$win_peak_info.="$win_peakNAM[$i]:chr$farm[0]:$win_peakstar[$i]-$win_peakend[$i]:$win_peaklod[$i];";
					if($win_peakNAM[$i]=~/NAM(\d+)/){
						push @NAM_id,$1;
					}
					push @win_peak_star,$win_peakstar[$i];
					push @win_peak_end,$win_peakend[$i];
					push @win_peak_lod,$win_peaklod[$i];
				}
			}
		}
		if($win_peak_info eq ""){
		}else{
			$QTL_num++;
			$QTL_mlm_info{$QTL_num}=$mlm_info;
			$QTL_win_peak_info{$QTL_num}=$win_peak_info;
			$QTL_peakP{$QTL_num}=$farm[2];
			$QTL_pos{$QTL_num}="Chr$farm[0]-$farm[1]";
			if($#win_peak_star==0){
				if($win_peak_star[0]>=200000){
					$star=$win_peak_star[0]-200000;
				}else{
					$star=0;
				}
				if($chr_length{$farm[0]}-$win_peak_end[0]>=200000){
					$end=$win_peak_end[0]+200000;
				}else{
					$end=$chr_length{$farm[0]}
				}
			}else{
				if((min @win_peak_lod)>=15){
					@new_win_peak_star = @win_peak_star;
					@new_win_peak_end = @win_peak_end;
				}else{
					foreach $i(0..$#win_peak_lod){
						if(($win_peak_lod[$i]/(max @win_peak_lod)) >= 0.5){
							push @new_win_peak_star,$win_peak_star[$i];
							push @new_win_peak_end,$win_peak_end[$i];
						}
					}
				}
				%duplicate_star = (); %duplicate_end = ();
				foreach (@new_win_peak_star){
					$duplicate_star{$_}++;
				}
				foreach (@new_win_peak_end){
					$duplicate_end{$_}++;
				}
				$num_star=0;$num_end=0;
				foreach $i(sort { $a <=> $b } @new_win_peak_star){
					if($duplicate_star{$i}>1){
						$star = $i;
						$num_star++;
						last;
					}
				}
				if($num_star==0){
					$star = min @new_win_peak_star;
				}
				foreach $i(sort { $b <=> $a } @new_win_peak_end){
					if($duplicate_end{$i}>1){
						$end = $i;
						$num_end++;
						last;
					}
				}
				if($num_end==0){
					$end = max @new_win_peak_end ;
				}
				if($star >= 200000){
					$star=$star-200000;
				}else{
					$star=0;
				}
				if($end+200000<=$chr_length{$farm[0]}){
					$end=$end+200000;
				}else{
					$end=$chr_length{$farm[0]};
				}
			}
			$geno_num=0;@farm_geno=();
			foreach $i($star..$end){
				if($position[$farm[0]][$i]=~/^LOC/ or $position[$farm[0]][$i]=~/^Os/ ){
					if(grep {$position[$farm[0]][$i] eq $_} @farm_geno){
					}else{
						$geno_num++;
						$QTL_farm_geno[$QTL_num][$geno_num]=$position[$farm[0]][$i];
						push @farm_geno,$position[$farm[0]][$i];
						push @geno,$position[$farm[0]][$i];
						@{$position[$farm[0]][$i]}=@NAM_id;
						$geno_id{$position[$farm[0]][$i]} ="Ture";
					}
				}
			}
			$qtl_genonum{$QTL_num}=$geno_num;
		}
	}
}
close IN9;
#open IN6,"SIFTtolerantScore.info" or die "can't open SIFTtolerantScore.info";
open IN6,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/Coding/SIFTtolerantScore_MSUvsRAP.txt" or die "can't open SIFTtolerantScore.info";
while (<IN6>){
	chomp();
	@sift=split/\s+/,$_;
	$key="$sift[0]	$sift[1]";
	$siftscore{$key}=1-$sift[5];
	if($sift[7]=~/STOP/ or $sift[7]=~/START/){
		$siftstop{$key}="Ture";
	}
	
}
close IN6;
foreach $chr(1..12){
	$file='/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/Coding/Chr'."$chr".'_combineSNPeffect_MSUvsRAP.txt';
	open IN,"$file" or die "can't open /data/home/mjchen/rice/NAM_denovo_assemble/mummer4/indel_blast/chr$chr-nam-Effect.snp";
	while (<IN>){
		chomp();
		@snp=split/\s+/,$_;
		if(defined($position[$chr][$snp[1]]) ){
			if($geno_id{$position[$chr][$snp[1]]} eq "Ture"){
				@genotype=();$mis=0;
				foreach $i(6..$#snp){
					if($i==12){
						if($snp[2] eq $snp[5]){
							push @genotype,0;
						}else{
							push @genotype,1;
						} 
						if($snp[$i] eq $snp[5]){
							push @genotype,0;
						}else{
							push @genotype,1;
						}
					}else{
						if($snp[$i] eq $snp[5]){
							push @genotype,0;
						}else{
							push @genotype,1;
						}
					}
				}
				foreach $i(@{$position[$chr][$snp[1]]}){
					if($genotype[$i-1] ==1){
						$mis++;
					}
				}
				
				if($mis>0){
					$key="$chr	$snp[1]	$snp[1]";
					$snp{$key}="Ture";
					$key="$chr	$snp[1]";
					if(defined($siftstop{$key})){
						$geneSTOP{$position[$chr][$snp[1]]}="T";
					}
				}
				$snpkey="$chr	$snp[1]";
				if(defined($siftscore{$snpkey})){
					foreach $i(5..$#snp){
						if($i==12){
							push @{"match$position[$chr][$snp[1]]"},0;
						}
						if($snp[$i] eq $snp[2]){
							push @{"match$position[$chr][$snp[1]]"},0;
						}else{
							push @{"match$position[$chr][$snp[1]]"},1;
						}
					}
					if($siftscore{$snpkey} <0.9){
						push @{"match_score$position[$chr][$snp[1]]"},1;
					}elsif($siftscore{$snpkey} <0.95 and $siftscore{$snpkey} >=0.9){
						push @{"match_score$position[$chr][$snp[1]]"},5;
					}elsif($siftscore{$snpkey} <0.98 and $siftscore{$snpkey} >=0.95){
						push @{"match_score$position[$chr][$snp[1]]"},15;
					}elsif($siftscore{$snpkey} <0.99 and $siftscore{$snpkey} >=0.98){
						push @{"match_score$position[$chr][$snp[1]]"},18;
					}elsif($siftscore{$snpkey}>= 0.99 and $siftscore{$snpkey}<1){
						push @{"match_score$position[$chr][$snp[1]]"},20;
					}elsif($siftscore{$snpkey}==1){
						if(defined $siftstop{$snpkey}){
							push @{"match_score$position[$chr][$snp[1]]"},22;
						}else{
							push @{"match_score$position[$chr][$snp[1]]"},20;
						}
					} 
					push @{"match_pos$position[$chr][$snp[1]]"},"SNP-$chr-$snp[1]";
				}
				#$match=0;
				#foreach $i(1..15){
				#	if(grep {$i eq $_}@{$position[$chr][$snp[1]]}){
				#		if($genotype[$i-1] ==1){
				#			$match++;
				#		}
				#	}else{
				#		if($genotype[$i-1] ==0){
				#			$match++;
				#		}
				#	}
				#}
				#$snpkey="$chr	$snp[1]";
				#if($match==15  and $snp[4] ne "*" and $siftscore{$snpkey}<=0.05){ 
				#	$aa=$position[$chr][$snp[1]];
				#	$match_score{$position[$chr][$snp[1]]}="ture";
				#	#print "aa	$snp[0]	$snp[1]	$position[$chr][$snp[1]]	@{$aa}	@genotype\n";
				#}	
			}
		}
	}
	close IN;
}


open IN6,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/Coding/SIFTtolerantScore_MSUvsRAP.txt" or die "can't open SIFTtolerantScore.info";
while (<IN6>){
	chomp();
	@sift=split/\s+/,$_;
	if($.==1){
		$key="$sift[0]	$sift[1]	$sift[2]";
		$siftscore{$key}=$sift[5];
		if($sift[6]=~/(Os\d+g\d+).\d+/){
			$geno=$1;
			if(defined($snp{$key})){
				
				if($sift[5]=~/\,/){
					@score_new=();
					@score_new=split/\,/,$sift[5];
					push @score,@score_new;
				}else{	
					push @score,$sift[5];
				}
			}
		}
	}else{
		$key="$sift[0]	$sift[1]	$sift[2]";
		if($sift[6]=~/(LOC_Os\d+g\d+).\d+/ or $sift[6]=~/(Os\d+g\d+).\d+/){
			if($sift[6]=~/(LOC_Os\d+g\d+).\d+/){
				$gene_nam=$1;
			}elsif($sift[6]=~/(Os\d+g\d+).\d+/){
				$gene_nam=$1;
			}
			if($geno eq $gene_nam){
				if(defined($snp{$key})){
					if($sift[5]=~/\,/){
						@score_new=();
						@score_new=split/\,/,$sift[5];
						push @score,@score_new;
					}else{	
						push @score,$sift[5];
					}
				}
			}else{
				if(!defined $score[0]){
					$sift{$geno}=0;
				}else{
					$sift{$geno}=1- (min @score);
				}
				$geno=$1;
				@score=();
				if(defined($snp{$key})){
					if($sift[5]=~/\,/){
						@score_new=();
						@score_new=split/\,/,$sift[5];
						push @score,@score_new;
					}else{	
						push @score,$sift[5];
					}
				}
			}
		}
	}
}  
if(!defined $score[0]){
	$sift{$geno}=0;	
}else{
	$sift{$geno}=1- (min @score);
}
@score=();
close IN6;
#非编码区snp变异打分
$hang=0;
open NoclodingSNP,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/Nocoding/combineRAPvsMSU_nocodingSNP.txt";
while (<NoclodingSNP>){
	chomp();
	$hang++;$mis=0;
	if($hang>1){
		@NoclodingSNP=split/\t+/,$_;
		foreach $i(@{$NoclodingSNP[2]}){
			if($NoclodingSNP[$i+3] ==1){
				$mis++;
				
			}
		}
		if($NoclodingSNP[11] ==0){
			push @{"match$NoclodingSNP[2]"},0;
			foreach $i(4..$#NoclodingSNP){
				push @{"match$NoclodingSNP[2]"},$NoclodingSNP[$i];
			}
		}else{
			push @{"match$NoclodingSNP[2]"},1;
			foreach $i(4..$#NoclodingSNP){
				if($NoclodingSNP[$i]==0){
					push @{"match$NoclodingSNP[2]"},1;
				}else{
					push @{"match$NoclodingSNP[2]"},0;
				}
			}
		}
		push @{"match_score$NoclodingSNP[2]"},$NoclodingSNP[3];
		push @{"match_pos$NoclodingSNP[2]"},"NoclodingSNP-$NoclodingSNP[0]-$NoclodingSNP[1]";
		if($NoclodingSNP[2] eq "LOC_Os01g66100"){
			$mmm="match$NoclodingSNP[2]";
			#print "NoclodingSNP	$NoclodingSNP[2]	$NoclodingSNP[3]	$NoclodingSNP[0]	$NoclodingSNP[1]	@{$mmm}\n";
		}
		if($mis > 0){
			push @{"Noclodingsnp$NoclodingSNP[2]"},$NoclodingSNP[3];
			if( grep {$NoclodingSNP[2] eq $_} @Noclodingsnp_gene){
			}else{
				push @Noclodingsnp_gene,$NoclodingSNP[2];
			}
		}
	}
}
close NoclodingSNP;
foreach $i(@Noclodingsnp_gene){
	$a="Noclodingsnp$i";
	$NoclodingSNPScore{$i}= max @{"Noclodingsnp$i"};
}
###indel变异打分信息 

foreach $chr(1..12){
	$file='/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/Coding/Chr'."$chr".'_combineINDELtolerantScore_MSUvsRAP.txt';
	open IN, "$file" or die "can't open /data/home/mjchen/rice/NAM_denovo_assemble/mummer4/indel_blast/chr$chr-nam-Effect.indel";
	while (<IN>){
		chomp();
		@indel=split/\s+/,$_;
		if(defined($position[$chr][$indel[1]]) or defined($position[$chr][$indel[2]])){
			if($geno_id{$position[$chr][$indel[1]]} eq "Ture" or $geno_id{$position[$chr][$indel[2]]} eq "Ture"){
				@genotype=();$mis=0;
				foreach $i(7..$#indel-1){
					if($i==13){
						if($indel[6]==0){
							push @genotype,0;
						}else{
							push @genotype,1;
						}
						if($indel[$i] eq $indel[6]){
							push @genotype,0;
						}else{
							push @genotype,1;
						}
					}else{
						if($indel[$i] eq $indel[6]){
							push @genotype,0;
						}else{
							push @genotype,1;
						}
					}
				}
				if($position[$chr][$indel[1]] eq $position[$chr][$indel[2]]){
					foreach $i(@{$position[$chr][$indel[1]]}){
						if($genotype[$i-1] ==1){
							$mis++;
						}
					}
				}elsif(defined($position[$chr][$indel[1]])){
					foreach $i(@{$position[$chr][$indel[1]]}){
						if($genotype[$i-1] ==1){
							$mis++;
						}
					}
				}elsif(defined($position[$chr][$indel[2]])){
					foreach $i(@{$position[$chr][$indel[2]]}){
						if($genotype[$i-1] ==1){
							$mis++;
						}
					}
				}
				#@{"indelScore$indel[5]"}=();
				
				if($mis>0){
					$key="Chr$chr	$indel[1]	$indel[2]";
					$indel{$key}="Ture"; 
					push @{"indelScore$indel[5]"},$indel[$#indel];
				}else{
					$key="Chr$chr	$indel[1]	$indel[2]"; 
				}
				
				 
				foreach $i(6..$#indel-1){
					if($i==13){
						push @{"match$position[$chr][$indel[1]]"},0;
					}
					push @{"match$position[$chr][$indel[1]]"},$indel[$i];
				}
				push @{"match_score$position[$chr][$indel[1]]"},$indel[$#indel];
				push @{"match_pos$position[$chr][$indel[1]]"},"INDEL-$chr-$indel[1]";
				
				#$match=0;
				#foreach $i(1..15){
				#	if(grep {$i eq $_}@{$position[$chr][$indel[1]]}){
				#		if($genotype[$i-1] ==1){
				#			$match++;
				#		}
				#	}else{
				#		if($genotype[$i-1] ==0){
				#			$match++;
				#		}
				#	}
				#}
				
				#if($match==15 and $indel[$#indel] >= 10){
				#	$aa="$position[$chr][$indel[1]]";
				#	$match_score{$position[$chr][$indel[1]]}="ture";
				#	#print "bb	$indel[0]	$indel[1]	$position[$chr][$indel[1]]	@{$aa}	@genotype\n";
				#}	
			}
		}
	}
	close IN;
}
#非编码区indel变异打分
$hang=0;
open NocodingIndel,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/Nocoding/combineRAPvsMSU_nocodingINDEL.txt";
while (<NocodingIndel>){
	chomp();
	$hang++;
	if($hang>1){
		my @NocodingIndel=split/\t+/,$_;
		
		@a=();@b=();@c=();
		if($NocodingIndel[3]=~/;/){
			@a=split/;/,$NocodingIndel[3];
			foreach $i(@a){
				if($i=~/(\S+):(\d+)/){
					push @b,$1;
					push @c,$2;
				}
			}
		}else{
			if($NocodingIndel[3]=~/(\S+):(\d+)/){
				push @b,$1;
				push @c,$2;
			}
		}
		foreach $j(0..$#b){
			$mis=0; 
			foreach $i(@{$b[$j]}){
				if($NocodingIndel[$i+3] ==1){
					$mis++;
				}
			}
			if($mis>0){
				push @{"Noclodingindel$b[$j]"},$c[$j];
				if( grep {$b[$j] eq $_} @Noclodingindel_gene){
				}else{
					push @Noclodingindel_gene,$b[$j];
				}
			}
			if($NocodingIndel[11] ==0){
				push @{"match$b[$j]"},0;
				foreach $i(4..$#NocodingIndel){
					push @{"match$b[$j]"},$NocodingIndel[$i];
				}
			}else{
				push @{"match$b[$j]"},1;
				foreach $i(4..$#NocodingIndel){
					if($NocodingIndel[$i]==0){
						push @{"match$b[$j]"},1;
					}else{
						push @{"match$b[$j]"},0;
					}
				}
			}
			
			
			
			
			push @{"match_score$b[$j]"},$c[$j];
			push @{"match_pos$b[$j]"},"NocodingIndel-$NocodingIndel[0]-$NocodingIndel[1]";
		}
	}
}
foreach $i(@Noclodingindel_gene){
	$a="Noclodingindel$i";
	$NoclodingIndelScore{$i}= max @{"Noclodingindel$i"};
}

foreach $chr(1..12){
	$file='/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/minigraph/NAM_Chr'."$chr".'combineRAPvsMSU_codingSV.txt';
	open IN8,"$file" or die "can't open NAM_Chr$chr-geneID.sv";
	while (<IN8>){
		chomp();
		@sv=split/\t+/,$_;
		if($.>1){
			if($_ eq "*"){
			}else{
				@sv_geno=();
				foreach $i(6..$#sv-1){
					if($i==12){
						if($sv[5] eq "0|0"){
							push @sv_geno,0;
						}else{
							push @sv_geno,1;
						}
						if($sv[5] eq $sv[$i]){
							push @sv_geno,0;
						}else{
							push @sv_geno,1;
						}
					}else{
						if($sv[5] eq $sv[$i]){
							push @sv_geno,0;
						}else{
							push @sv_geno,1;
						}
					}
				}
				foreach $i(@geno){
					if($sv[0]=~/$i/){
						$miss=0;
						foreach $j(@{$i}){
							if($sv_geno[$j-1] == 1){
								$miss++;
							}
						}
						if($miss>0){
							push @{"svScore$sv[0]"},$sv[$#sv];
							$sv_score{$i}="1";
						}else{
							$sv_score{$i}="*";
						}
					}
					 
				}
				
				foreach $i(5..$#sv-1){
					if($i==12){
						push @{"match$sv[0]"},0;
					}
					if($sv[$i] eq "0|0"){
						push @{"match$sv[0]"},0;
					}else{
						push @{"match$sv[0]"},1;
					}
				}
				push @{"match_score$sv[0]"},$sv[$#sv];
				push @{"match_pos$sv[0]"},"SV-$sv[1]-$sv[2]";
				if($sv[0] eq "LOC_Os01g66100"){
					$mmm="match$sv[0]";
					#print "SV	$sv[0]	$sv[$#sv]	$sv[1]	$sv[2]	@{$mmm}\n";
				}
				
				#$match=0;
				#foreach $i(1..15){
				#	if(grep {$i eq $_}@{$sv[0]}){
				#		if($sv_geno[$i-1] ==1){
				#			$match++;
				#		}
				#	}else{
				#		if($sv_geno[$i-1] ==0){
				#			$match++;
				#		}
				#	}
				#}
				#if($match==15){
				#	$match_score{$sv[0]}="ture";
					#print "cc	$sv[1]	$sv[2]	$sv[0]	@{$sv[0]}	@sv_geno\n";
				#}
			}
		}
	}
	close IN8;
}

#非编码区sv变异打分
$hang=0;
open NocodingSV,"/data/home/mjchen/rice/NAM_denovo_assemble/Nested_association_mapping/illumina_data/NAM-all-info/QTLtoPearson/Nocoding/combineRAPvsMSU_nocodingSV.txt";
while (<NocodingSV>){
	chomp();
	$hang++;
	if($hang>1){
		my @NocodingSV=split/\t+/,$_;
		@a=();@b=();@c=();
		if($NocodingSV[3]=~/;/){
			@a=split/;/,$NocodingSV[3];
			foreach $i(@a){
				if($i=~/(\S+):(\d+)/){
					push @b,$1;
					push @c,$2;
				}
			}
		}else{
			if($NocodingSV[3]=~/(\S+):(\d+)/){
				push @b,$1;
				push @c,$2;
			}
		}
		foreach $j(0..$#b){
			$mis=0; 
			foreach $i(@{$b[$j]}){
				if($NocodingSV[$i+3] ==1){
					$mis++;
				}
			}
			if($mis>0){
				push @{"NoclodingSV$b[$j]"},$c[$j];
				if( grep {$b[$j] eq $_} @NoclodingSV_gene){
				}else{
					push @NoclodingSV_gene,$b[$j];
				}
			}
			 
			if($NocodingSV[11] ==0){
				push @{"match$b[$j]"},0;
				foreach $i(4..$#NocodingSV){
					push @{"match$b[$j]"},$NocodingSV[$i];
				}
			}else{
				push @{"match$b[$j]"},1;
				foreach $i(4..$#NocodingSV){
					if($NocodingSV[$i]==0){
						push @{"match$b[$j]"},1;
					}else{
						push @{"match$b[$j]"},0;
					}
				}
			}
			push @{"match_score$b[$j]"},$c[$j];
			push @{"match_pos$b[$j]"},"NocodingSV-$NocodingSV[0]-$NocodingSV[1]";
			if($b[$j] eq "LOC_Os01g66100"){
				$mmm="match$b[$j]";
				#print "NoclodingSV	$b[$j]	$c[$j]	$NocodingSV[0]	$NocodingSV[1]	@{$mmm}\n";
			}
		}
	}
}
foreach $i(@NoclodingSV_gene){
	$a="NoclodingSV$i";
	$NoclodingSVScore{$i}= max @{"NoclodingSV$i"};
}



@ranseq_name=("root_seeding_stage","shoot_seeding_stage","leaf_sheath_seeding_stage","pistils","anthers","Leaves_seeding_stage","leaves_tillering_stage","leaves_flowering_stage","young_panicle","panicle_flowering_stage","panicle_filling_stage","tiller_buds","embryos","developing_seeds","seeds_after_ageing","seed_Germinating_stage","shoot_apical_meristem","lamina_joint");

foreach $j(0..$#ranseq_name){
	if($ranseq_name[$j] eq $ARGV[1]){
		$main_RNA=$j;
	}
}
open OUT,">RAP_$ARGV[0].FarmCPUpeak_info";
print OUT "QTL	FarmCPU_pos	FarmCPU_-logP	MLM_INFO(Chr-pos-(-log10)	Winqtlcart_peak_info	Geno	Geno_pos	LOD	transposon	House_keeping	Importance	AT	RAPDB_annotation	SIFT_Score	INDEL_Score	SV_Score	Match_Score	root_seeding_stage	shoot_seeding_stage	leaf_sheath_seeding_stage	pistils	anthers	Leaves_seeding_stage	leaves_tillering_stage	leaves_flowering_stage	young_panicle	panicle_flowering_stage	panicle_filling_stage	tiller_buds	embryos	developing_seeds	seeds_after_ageing	seed_Germinating_stage	shoot_apical_meristem	lamina_joint\n";
#print "Gene_id	match_time	type	P1	P2	P3	P4	P5	P6	P7	P8	P9	P10	P11	P12	P13	P14	P16	P17	winQtlcart\n";
foreach $m(1..$QTL_num){			
	foreach $n(1..$qtl_genonum{$m}){
		$QTL=$m	;
		$i=$QTL_farm_geno[$m][$n];
		#if($annotation{$i} ne "transposon" and  $annotation{$i} ne "retrotransposon" and $house_keeping{$i} ne "House_keeping"){
			#if($sift{$i} >= 0.95 or (max(@{"indelScore$i"})) >= 5 or (max(@{"svScore$i"})) >= 5){
			#if($sift{$i} >= 0 or (max(@{"indelScore$i"})) >= 0 or (max(@{"svScore$i"})) >= 0){
				@geno_rna=split/\s+/,$rna_geno{$i};
				$highsampleNum=0;
				$RNAselect="False";
				if($geno_rna[$main_RNA+1] > 0 and $geno_rna[$main_RNA+1]<1){
					foreach $j(1..$#geno_rna){
						if($j==$main_RNA+1){
						}else{
							if($geno_rna[$j] >= 50*$geno_rna[$main_RNA+1]){
								$highsampleNum++;
							}
						}
					}
					if($highsampleNum <=5){
						$RNAselect="Ture";
					}else{
						$RNAselect="False";
					}
				}elsif($geno_rna[$main_RNA+1]>=1){
					foreach $j(1..$#geno_rna){
						if($j==$main_RNA+1){
						}else{
							if($geno_rna[$j] >= 20*$geno_rna[$main_RNA+1]){
								$highsampleNum++;
							}
						}
					}
					if($highsampleNum <=5){
						$RNAselect="Ture";
					}else{
						$RNAselect="False";
					}
				}
				#if($RNAselect eq "Ture"){
					print OUT "QTL$QTL	$QTL_pos{$m}	$QTL_peakP{$m}	$QTL_mlm_info{$m}	$QTL_win_peak_info{$m}	$i";
					$lod_name="LOD$i";
					print OUT "	$pos{$i}	@{$lod_name}";
					
					if (defined ($annotation{$i})){
						print OUT "	$annotation{$i}";
					}else{
						print OUT "	*";
					}
					if (defined ($house_keeping{$i})){
						print OUT "	$house_keeping{$i}";
					}else{
						print OUT "	*";
					}
					if (defined ($Important{$i})){
						print OUT "	$Important{$i}";
					}else{
						print OUT "	*";
					}
					if (defined ($AT{$i})){
						print OUT "	$AT{$i}";
					}else{
						print OUT "	*";
					}
					if (defined ($msu_rap_ann{$i})){
						print OUT "	$msu_rap_ann{$i}";
					}else{
						print OUT "	*";
					}
					if (defined ($sift{$i})){
						if($sift{$i} <0.9){
							$sift_snp_score=1;
						}elsif($sift{$i} <0.95 and $sift{$i} >=0.9){
							$sift_snp_score=5;
						}elsif($sift{$i} <0.98 and $sift{$i} >=0.95){
							$sift_snp_score=15;
						}elsif($sift{$i} < 0.99 and $sift{$i} >=0.98){
							$sift_snp_score=18;
						}elsif($sift{$i} >= 0.99){
							$sift_snp_score=20;
						}
						if(defined ($geneSTOP{$i}) ){
							$sift_snp_score=25;
						}
						if(defined($NoclodingSNPScore{$i})){
							if($sift_snp_score > $NoclodingSNPScore{$i}){
								print OUT "	$sift_snp_score";
							}else{
								print OUT "	$NoclodingSNPScore{$i}";
							}
						}else{
							print OUT "	$sift_snp_score";
						}
					}else{
						if(defined ($NoclodingSNPScore{$i})){
							print OUT "	$NoclodingSNPScore{$i}"
						}else{
							print OUT "	*";
						}	
					}
					if (defined (${"indelScore$i"}[0])){
						$indel_score=max(@{"indelScore$i"});
						
						if(defined ($NoclodingIndelScore{$i})){
							if($indel_score>$NoclodingIndelScore{$i}){
								print OUT "	$indel_score";
							}else{
								print OUT "	$NoclodingIndelScore{$i}";
							}
						}else{
							print OUT "	$indel_score";
						}
					}else{
						if(defined ($NoclodingIndelScore{$i})){
							print OUT "	$NoclodingIndelScore{$i}"
						}else{
							print OUT "	*";
						}
					}
					if (defined (${"svScore$i"}[0])){
						$sv_score=max(@{"svScore$i"});
						if(defined ($NoclodingSVScore{$i})){
							if($sv_score > $NoclodingSVScore{$i}){
								print OUT "	$sv_score";
							}else{
								print OUT "	$NoclodingSVScore{$i}";
							}
						}else{
							print OUT "	$sv_score";
						}
					}else{
						if(defined ($NoclodingSVScore{$i})){
							print OUT "	$NoclodingSVScore{$i}"
						}else{
							print OUT "	*";
						}
					}
					@ture_geno=(); 
					$match_s=0;
					$a="match$i";
					$b="match_score$i";
					$c="match_pos$i";
					#print "aa	$i	$#{$a}	$#{$b}	@{$b}	@{$a}\n";
					if(max @{$b} >=20){
						$match_num=0;
						foreach $m_score(0..$#{$b}){
							if(${$b}[$m_score] >= 22){
								$match_num++;
								$aaa=$m_score*16;
								$bbb=($m_score+1)*16-1;
								#print "$i	$match_num	${$c}[$m_score]";
								foreach $m_gene($m_score*16..($m_score+1)*16-1){
									#print "	${$a}[$m_gene]";
									if(${$a}[$m_gene] eq 1){
										$order=$m_gene-$m_score*16+1;
										if(grep {$order eq $_}@ture_geno ){
										}else{
											push @ture_geno,$order;
										}
									}
								}
								#print "	@{$i}\n";
							}
						}
						
						foreach $e(1..16){
							if(grep {$e eq $_}@ture_geno){
								$geno_base[$e-1]=1;
							}else{
								$geno_base[$e-1]=0;
							}
						}
						@ture_genobase=();
						foreach $e(1..15){
							if($geno_base[$e] eq $geno_base[0]){
							}else{
								push @ture_genobase,$e;
							}
						}
						$q=0;
						
						if($#ture_genobase eq $#{$i}){
							foreach $ture_g(@ture_genobase){
								if(grep {$ture_g eq $_} @{$i}){
									$q++;
								}
							}
							if ($q == $#{$i}+1){
								$match_s=5;
							}
						}
						
					}elsif(max @{$b} >=10 and max @{$b} <=20){ 
						$match_num=0;
						foreach $m_score(0..$#{$b}){
							if( ${$b}[$m_score] eq max @{$b}){
								$match_num++;
								@ture_geno=();
								$aaa=$m_score*16;
								$bbb=($m_score+1)*16-1;
								#print "$i	$match_num	${$c}[$m_score]";
								foreach $m_gene($m_score*16..($m_score+1)*16-1){
									#print "	${$a}[$m_gene]";
									if(${$a}[$m_gene] eq 1){
										$order=$m_gene-$m_score*16+1;
										if(grep {$order eq $_}@ture_geno ){
										}else{
											push @ture_geno,$order;
										}
									}
								}
								#print "	@{$i}\n";
							}
						}
						foreach $e(1..16){
							if(grep {$e eq $_}@ture_geno){
								$geno_base[$e-1]=1;
							}else{
								$geno_base[$e-1]=0;
							}
						}
						@ture_genobase=();
						foreach $e(1..15){
							if($geno_base[$e] eq $geno_base[0]){
							}else{
								push @ture_genobase,$e;
							}
						}
						$q=0;
						if($#ture_genobase eq $#{$i}){
							foreach $ture_g(@ture_genobase){
								if(grep {$ture_g eq $_} @{$i}){
									$q++;
								}
							}
							if ($q == $#{$i}+1 ){
								if($match_num==1){
									$match_s=5;
								}else{
									$match_s="abnormal";
								}
							}
						}
					}
					
					if($match_s==5){
						print OUT "	5";
					}elsif($match_s==0){
						print OUT "	*";
					}elsif($match_s eq "abnormal"){
						print OUT "	abnormal";
					}
						
					if(defined($rna_geno{$i})){
						print OUT "$rna_geno{$i}";
					}else{
						print OUT "	*";
					}
					
					
					print OUT "\n";
				#}
			#}
		#}
	}
}
close OUT;
close IN9;
			
		
		
