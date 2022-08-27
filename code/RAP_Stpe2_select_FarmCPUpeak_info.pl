#perl Stpe2_select_FarmCPUpeak_info.pl SH_Plant_height shoot_seeding_stage
use List::Util qw/min sum max/;
open IN1,"key_word.txt" or die "can't open key_word.txt";
while (<IN1>){
	chomp();
	if($_=~/(\S+):(.*);(.*)?/){
		$trait=$1;
		$key=$2;
		$key2=$3;
		if($ARGV[1]=~/$trait/i){
			@key_words=split/\,/,$key;
			@key_words_masked=split/\,/,$key2;
			
		}
	}
}
close IN1;
$tissue{root_seeding_stage}="";	
$tissue{shoot_seeding_stage}="Leaf_length,Leaf_width,Plant_height,Culm_length,Panicle_enclosure,Panicle_length";
$tissue{leaf_sheath_seeding_stage}="";
$tissue{pistils}="Recombination_times_per_sample";
$tissue{anthers}="Recombination_times_per_sample";
$tissue{Leaves_seeding_stage}="Leaf_length,Leaf_width";
$tissue{leaves_tillering_stage}="";
$tissue{leaves_flowering_stage}="Leaf_angle,Leaf_length,Leaf_width,Heading_data,Plant_height,Culm_length,Panicle_enclosure,Protein,Seed_length,Seed_width";
$tissue{young_panicle}="Hull_color,Awn_length,Panicle_length,Yield,Seed_length,Seed_width";
$tissue{panicle_filling_stage}="Panicle_length,Panicle_enclosure";
$tissue{tiller_buds}="Tiller_number,Mutiple_panicles_per_tiller,Yield";
$tissue{embryos}="";
$tissue{developing_seeds}="Protein";
$tissue{seeds_after_ageing}="";
$tissue{seed_Germinating_stage}="48H,60H,72H";
$tissue{shoot_apical_meristem}="Leaf_length,Leaf_width,Plant_height,Culm_length,Panicle_enclosure,Panicle_length";
$tissue{lamina_join}="Leaf_angle";

@ranseq_name=("root_seeding_stage","shoot_seeding_stage","leaf_sheath_seeding_stage","pistils","anthers","Leaves_seeding_stage","leaves_tillering_stage","leaves_flowering_stage","young_panicle","panicle_flowering_stage","panicle_filling_stage","tiller_buds","embryos","developing_seeds","seeds_after_ageing","seed_Germinating_stage","shoot_apical_meristem","lamina_joint");
foreach $j(0..$#ranseq_name){
	if($tissue{$ranseq_name[$j]}=~ /$ARGV[1]/){
		push @main_RNA,$j+17;
	}
}
@farmScore=();
@mlmScore=();
@WinqtlcartScore=();
open IN2,"RAP_$ARGV[0].FarmCPUpeak_info" or die "can't open $ARGV[0].FarmCPUpeak_info";
while (<IN2>){
	chomp(); 
	if($.>1){
		$line=$_;
		@line =split/\t+/,$line;
		if(grep {$line[2] eq $_}@farmScore){
		}else{
			push @farmScore,$line[2];
		}
		@wininfo=split/\;/,$line[4];
		foreach $i(0..$#wininfo){
			if($wininfo[$i]=~/NAM\d+\:chr\d+:\d+\-\d+\:(\S+)/){
				$score=$1;
				if(grep {$score eq $_} @WinqtlcartScore){
				}else{
					push @WinqtlcartScore,$score;
				}
			}
		}
	}	
}
close IN2;		
@sort_farmScore = sort {$a <=> $b}  @farmScore; 
@sort_WinqtlcartScore = sort  {$a <=> $b} @WinqtlcartScore;
if($#sort_farmScore>=4){
	$rank_farm=2;
	$group_farm=5;
}elsif($#sort_farmScore<=3){
	
	$rank_farm=10/($#sort_farmScore+1);
	$group_farm=$#sort_farmScore+1;
}
if($#sort_WinqtlcartScore>=4){
	$rank_Win=2;
	$group_Win=5;
}elsif($#sort_WinqtlcartScore<=3){
	$rank_Win=10/($#sort_WinqtlcartScore+1);
	$group_Win=$#sort_WinqtlcartScore+1;
}


open OUT,">RAP_Step2_$ARGV[0].FarmCPUpeak_info";
print OUT "QTL	FarmCPU_pos	FarmCPU_-logP	MLM_INFO(Chr-pos-(-log10))	Winqtlcart_peak_info	Geno	Geno_pos	LOD	transposon	House_keeping	Importance	AT	RAPDB_annotation	SIFT_Score	INDEL_Score	SV_Score	Score_peak	Score_TE	Score_annotate	Score_expression	Score_match	Score	root_seeding_stage	shoot_seeding_stage	leaf_sheath_seeding_stage	pistils	anthers	Leaves_seeding_stage	leaves_tillering_stage	leaves_flowering_stage	young_panicle	panicle_flowering_stage	panicle_filling_stage	tiller_buds	embryos	developing_seeds	seeds_after_ageing	seed_Germinating_stage	shoot_apical_meristem	lamina_joint\n";
open IN2,"RAP_$ARGV[0].FarmCPUpeak_info" or die "can't open $ARGV[0].FarmCPUpeak_info";
while (<IN2>){
	chomp();
	if($.>1){
		$line=$_;
		@line =split/\t+/,$line;
		#farmCPU打分
		foreach $i(0..$#sort_farmScore){
			if($line[2] eq $sort_farmScore[$i]){
				$group_seq=int (($i)/(int(($#sort_farmScore+1)/$group_farm)));
				if($group_seq<=$group_farm){
					$score1_1=($group_seq+1)*$rank_farm;
				}elsif($group_seq>$group_farm){
					$score1_1=10;
				}
				last;
			}
		}
		#mlm打分
		if($line[3]=~/Chr\d+\-\d+\-(\S+)/){
			if($1<10){
				$score1_2=0;
			}elsif($1>=10 and $1<20){
				$score1_2=5;
			}elsif($1>=20){
				$score1_2=10;
			}
		}
		#Winqtlcart打分
		@wininfo=();
		@wininfo=split/\;/,$line[4];@s_win=();
		foreach $i(0..$#wininfo){
			if($wininfo[$i]=~/NAM\d+\:chr\d+:\d+\-\d+\:(\S+)/){
				$score_win=$1;
				foreach $j(0..$#sort_WinqtlcartScore){
					if($score_win eq $sort_WinqtlcartScore[$j]){
						$group_seq=int (($j)/(int(($#sort_WinqtlcartScore+1)/$group_Win)));
						if($group_seq<=$group_Win){
							$s_win[$i]=($group_seq+1)*$rank_Win;
						}elsif($group_seq>$group_Win){
							$s_win[$i]=10;
						}
						last;
						
					}
				}
			}
		}
		$score1_3 = ((sum @s_win)/($#s_win+1));
		$score1=($score1_1+$score1_2+$score1_3);
		#转座子等打分
		$trans=();$house=();
		if($line[8] eq "transposon" or $line[8] eq "retrotransposon"){
			$trans=0;
		}else{
			$trans=1;
		}
		if($line[9] eq "House_keeping"){
			$house=0;
		} else{
			$house=1;
		}
		if($line[10] eq "Important"){
			$Im=5;
		}else{
			$Im=0;
		}
		$score2=$Im;
		#拟南芥注释打分
		$num_AT=0;
		$num_rice=0;
		@AT=split/=>/,$line[11];
		#print "$line[5]	";
		foreach $i(@key_words){
			$AT_info=();
			$AT_info="$AT[1]$AT[3]$AT[4]";
			if($AT_info=~/$i/){
				$num_AT++;
			}
			if($line[12]=~/$i/){
				$num_rice++;
			}
		}
		foreach $i(@key_words_masked){
			$AT_info=();
			$AT_info="$AT[1]$AT[3]$AT[4]";
			if($AT_info=~/$i/){
				$num_AT--;
			}
			if($line[12]=~/$i/){
				$num_rice--;
			}
		}
		if($num_AT >=2 and $num_rice >=2){
			$score3=25;
		}elsif($num_AT >=3 and $num_rice ==1){
			$score3=23;
		}elsif($num_rice >=3 and $num_AT ==1){
			$score3=23;
		}elsif($num_AT >=3 and $num_rice ==0){
			$score3=18;
		}elsif($num_rice >=3 and $num_AT ==0){
			$score3=18;		
		}elsif($num_AT >=2 and $num_rice ==1){
			$score3=22;
		}elsif($num_rice >=2 and $num_AT ==1){
			$score3=22;
		}elsif($num_AT >=2 and $num_rice ==0){
			$score3=17;
		}elsif($num_rice >=2 and $num_AT ==0){
			$score3=17;	
		}elsif($num_AT ==1 and $num_rice >=1){
			$score3=20;
		}elsif($num_AT ==1 and $num_rice >=0){
			$score3=10;
		}elsif($num_AT ==0 and $num_rice >=1){
			$score3=10;
		}elsif($num_AT ==0 and $num_rice ==0){
			$score3=0;
		} 
		#表达量打分
		@select_main_RNA=0;
		foreach $m(@main_RNA){
			push @select_main_RNA,$line[$m];
		}
		foreach $m(@main_RNA){
			if($line[$m] == max @select_main_RNA){
				$main_RNA=$m;
			}
		}
			
		$rna=0;@rna=();@sort_rna=();$num=0;
		if($line[$main_RNA]>=0.05){
			$score4=10;
		}else{
			$score4=0;
		}
		foreach $i(17..$#line){
			if($i == $main_RNA){
			}else{
				push @rna,$line[$i];
			}
		}
		@sort_rna = sort {$b <=> $a}  @rna; 
		foreach (@sort_rna){
			if(($_/10)>$line[$main_RNA]){
				$num++;
			}
		}
		if($num>=5){
			$score4=5;
		}else{
			shift (@sort_rna);shift (@sort_rna);
			$average_rna=((sum @sort_rna)/($#sort_rna+1));
			if($line[$main_RNA]>0.05){
				if($average_rna==0){
					$average_rna=0.01;
				}
				if($line[$main_RNA]/$average_rna>=10){
					$score4=22;
				}elsif($line[$main_RNA]/$average_rna>=5){
					$score4=20;
				}elsif($line[$main_RNA]/$average_rna>=2){
					$score4=18;
				}elsif($line[$main_RNA]/$average_rna>=1){
					$score4=15;
				}else{
					$score4=10;
				}
			}
		}
		#if($trans==0 or $house==0 ){去掉housekeeping gene打分为0的规则
		if($trans==0){	
			$score=0;
		}else{
			@var=();
			if($line[13] ne "*"){
				push @var,$line[13];
			}
			if($line[14] ne "*"){
				push @var,$line[14];
			}
			if($line[15] ne "*"){
				push @var,$line[15];
			}
			if($#var>=0){
				$max=max @var;
			}else{
				$max=0;
			}
			if($line[16] ne "*"){	
				$score=$score1+$score2+$score3+$score4+$max+$line[16];
			}else{
				$score=$score1+$score2+$score3+$score4+$max;
			}
		}
		foreach (0..15){
			print OUT "$line[$_]	";
		}
		print OUT "$score1	$score2	$score3	$score4	$line[16]	$score";
		foreach (17..$#line){
			print OUT "	$line[$_]";
		}
		print OUT "\n";
	}
}
close IN2;
close OUT;		