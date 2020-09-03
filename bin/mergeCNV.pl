use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use List::Util qw(uniq);
use FindBin qw/$Bin/;
use Data::Dumper;


my ($sampleName, $infile, $bed, $annot_file, $strand_file, $gene_cutoff, $rpt_gene_file, $outdir);

GetOptions(
	"n:s" => \$sampleName,                             # need
	"in:s" => \$infile,                                # need
	"bed:s" => \$bed,                                  # need
	"annot:s" => \$annot_file,                         # need
	"strand:s" => \$strand_file,                       # need
	"pct:f" => \$gene_cutoff,                          # default: 0.7
	"rpt:s" => \$rpt_gene_file,                        # need
	"o:s" => \$outdir,                                 # default: dirname(infile)
	) or die "unknown args\n";

# check args
if (not defined $sampleName || not defined $infile || not defined $bed || not defined $annot_file || not defined $strand_file){
	die "please check your args\n";
}

# default value
if (not defined $gene_cutoff){
	$gene_cutoff = 0.7;
}

if (not defined $outdir){
	$outdir = dirname($infile);
}


# check if file exists
if (!-e $annot_file || !-e $bed || !-e $strand_file || !-e $infile){
	die "can not find $annot_file file or $bed file or $strand_file or $infile\n";
}

# ttranscript info
my $base_dir = dirname($Bin);
my $main_NM_file = "$base_dir/public_db/MainNM_yrt_20191111.txt";

# 结果文件
my $exon_level_res = "$outdir/$sampleName\.exon.CNV.xls";
my $gene_level_res = "$outdir/$sampleName\.gene.CNV.xls";

open O1, ">$exon_level_res" or die;
print O1 "sample\tgene\tchr\tstart\tend\tcopyNumber\tcnv_state\ttranscript\texon\n";

open O2, ">$gene_level_res" or die;
print O2 "sample\tgene\tchr\tstart\tend\tcopyNumber\tcnv_state\ttranscript\n";


# 基因+/-链信息
my %gene_strand;
open IN, "$strand_file" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$gene_strand{$arr[3]} = $arr[-4];
}
close IN;

# 基因坐标等信息
my (%gene_targetNum, %gene_pos, %gene_chr, @gene);

open BED, "$bed" or die;
while (<BED>){
	chomp;
	my @arr = split /\t/;
	next if ($arr[-1] eq "-");
	next if ($arr[-1] eq "HIST1H1C"); # no use gene
	$gene_chr{$arr[3]} = $arr[0]; # 基因在哪条染色体
	$gene_targetNum{$arr[3]}++; # 几个捕获区域
	push @{$gene_pos{$arr[3]}},$arr[1]; # start pos
	push @{$gene_pos{$arr[3]}},$arr[2]; # end pos
	push @gene,$arr[3]; # gene name
}
close BED;

# 外显子水平基因列表
my @exon_level_gene; # BRCA1 BRCA2 and other genes
open IN, "$rpt_gene_file" or die;
while (<IN>){
	chomp;
	my @arr = split /\,/, $_;
	push @exon_level_gene, $arr[0];
}

my %exon_level_gene;
for my $g (@exon_level_gene){
	$exon_level_gene{$g} = 1;
}
undef @exon_level_gene;


# 去重
my @gene_by_chr_order;
my %gene;
for my $g (@gene){
	if (!exists $gene{$g}){
		push @gene_by_chr_order, $g;
		$gene{$g} = 1;
	}else{
		next;
	}
}

#my @gene_by_chr_order = uniq(@gene);
undef @gene;

# 位置信息
my %gene_pos_sp;
for my $g (@gene_by_chr_order){
	my @pos = @{$gene_pos{$g}};
	$gene_pos_sp{$g}{"start"} = $pos[0];
	$gene_pos_sp{$g}{"end"} = $pos[-1];
}


######################################## CNV分析 ####################################

# for exon level
my %exon_level_gene_target_idx;
open IN, "$infile" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $target = "$arr[1]\t$arr[2]\t$arr[3]"; # chr/start/end
	if (exists $exon_level_gene{$arr[4]}){
		# gene=>{target=>"ratio/sampleQC/targetQC/CNV_state"}
		$exon_level_gene_target_idx{$arr[4]}{$target} = "$arr[9]\t$arr[12]\t$arr[11]\t$arr[-4]"; # ratio/sampleQC/targetQC/CNV_state [ratio/targetQC/CNV_state可能为NA]
	}
}
close IN;

# for gene level
my $target_idx;
my %gene_level_gene_target_idx;
open IN, "$infile" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	if (not exists $gene_level_gene_target_idx{$arr[4]}){
		$target_idx = 1;
		$gene_level_gene_target_idx{$arr[4]}{$target_idx} = "$arr[9]\t$arr[12]\t$arr[11]\t$arr[-4]";
	}else{
		$target_idx += 1;
		$gene_level_gene_target_idx{$arr[4]}{$target_idx} = "$arr[9]\t$arr[12]\t$arr[11]\t$arr[-4]";
	}
}
close IN;
undef $target_idx;

# 捕获区间对应的外显子编号
my %target2exon;
open IN, "$annot_file" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $target = "$arr[0]\t$arr[1]\t$arr[2]"; # chr/start/end
	$target2exon{$arr[3]}{$target} = $arr[-2]; # gene=>{target=>exon}
}
close IN;

# 外显子对应的捕获区域
my %exon2region;
open IN, "$annot_file" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $target = "$arr[0]\t$arr[1]\t$arr[2]"; # chr/start/end
	$exon2region{$arr[3]}{$arr[-2]} = $target; # gene=>{exon=>target}
}
close IN;


# 转录本号信息
my $NM_href = &get_main_NM($main_NM_file);

# 循环每一个基因
for my $gene (@gene_by_chr_order){
	my $NM; # 转录本号
	if (exists $NM_href->{$gene}){
		$NM = $NM_href->{$gene};
	}else{
		print "can not find NM for gene, will be skipped\n";
		next;
	}

	next if ($gene eq "SCG5" || $gene eq "GREM1"); # 跳过这2个基因,因为GREM1启动子区包括这2基因的绝大部分
	my $exon_num = $gene_targetNum{$gene}; # 外显子数目
	my $chr = $gene_chr{$gene}; # 染色体
	my $start = $gene_pos_sp{$gene}->{"start"}; # 起始坐标
	my $end = $gene_pos_sp{$gene}->{"end"}; # 终止坐标
	my $strand = $gene_strand{$gene}; # +/-链
	my $target = "$chr\t$start\t$end";

	if (exists $exon_level_gene{$gene}){
		# 外显子水平CNV
		print "#exon-level# $gene info:\n";
		my $gene_info_href = $exon_level_gene_target_idx{$gene}; # ratio/cnv_state等信息
		my @target = keys %{$gene_info_href};
		
		my %cnv_state; # 该基因每个外显子的cnv_state信息
		for my $t (@target){
			# 循环每个外显子
			if (exists $target2exon{$gene}{$t}){
				# skip intron region
				my $val = $gene_info_href->{$t};
				my $cnv_state = (split /\t/, $val)[-1];
				my $e = $target2exon{$gene}{$t}; # get exon number from annot file, will skip intron region
				$e =~ s/^exon//;
				$cnv_state{$e} = $cnv_state; # exon=>Del/Amp/Normal/...
			}
		}

		# 打印每个外显子cnv_state
		print(Dumper(\%cnv_state));

		my @type = qw/Amp Del/;
		for my $t (@type){
			# 分别输出Amp/Del类型CNV信息
			my @exon_sort;
			for my $e (sort {$a <=> $b} keys %cnv_state){ # 外显子由1->N排序
				my $state = $cnv_state{$e};
				if ($state eq $t){
					push @exon_sort, $e;
				}
			}

			if (scalar(@exon_sort) == 0){
				# 不存在Del/Amp
				push @exon_sort, "NA"; # for compatibility with &get_cnv_idx func for gene level cnv
			}

			my $ok_idx_aref = \@exon_sort;
			print "$gene $t index is: @{$ok_idx_aref}\n";

			# 合并相邻exon
			my $seg_aref = &merge_single_exon_into_large_seg($ok_idx_aref);
			print "$gene $t merged single exon seg is: @{$seg_aref}\n";

			next if (@{$seg_aref} == 1 and $seg_aref->[0] eq "NA"); # no idx

			################################## 输出结果 ##################################
			for my $seg (@{$seg_aref}){
				my @exon; # 外显子号
				my $ratio_sum; # 计算平均ratio
				my $n; # 外显子数目
				my @pos; # 位置信息

				if ($seg =~ /\_/){
					my @arr = split /\_/, $seg;
					for my $i (@arr){
						$n += 1;
						my $exon = "exon".$i;
						my $target = $exon2region{$gene}{$exon};
						my $val = $exon_level_gene_target_idx{$gene}{$target}; # 外显子ratio信息 <ratio/sampleQC/targetQC/CNV_state>
						my @val = split "\t", $val;
						my $ratio = $val[0];
						$ratio_sum += $ratio;
						
						my @region = split /\t/, $exon2region{$gene}{$exon}; # chr/start/end
						my $sp = $region[1];
						my $ep = $region[2];
						$exon =~ s/^exon//;
						
						push @pos, $sp;
						push @pos, $ep;
						push @exon, $exon;
					}
				}else{
					$n = 1;
					my $exon = "exon".$seg;
					my $target = $exon2region{$gene}{$exon};
					my $val = $exon_level_gene_target_idx{$gene}{$target};
					my @val = split "\t", $val;
					my $ratio = $val[0];
					$ratio_sum = $ratio;
					
					
					my @region = split /\t/, $exon2region{$gene}{$exon}; # chr/start/end
					my $sp = $region[1];
					my $ep = $region[2];
					$exon =~ s/^exon//;

					push @pos,$sp;
					push @pos,$ep;
					push @exon,$exon;
				}

				# 输出这个seg信息
				my $avg_ratio = sprintf "%.2f", $ratio_sum/$n;
				my $CN = sprintf "%.2f", 2*$avg_ratio; # 拷贝数
				my $CNType = &define_CN_state($CN); # 杂合/纯合

				my $exon;
				if (@exon == 1){
					$exon = "Exon".$exon[0];
				}else{
					$exon = "Exon".$exon[0]."-".$exon[-1];
				}

				my @pos_sort = sort {$a <=> $b} @pos;
				my $start_exon = $pos_sort[0];
				my $end_exon = $pos_sort[-1];

				print O1 "$sampleName\t$gene\t$chr\t$start_exon\t$end_exon\t$CN\t$CNType\t$NM\t$exon\n";
			}
		}
	}else{
		# gene-level CNV
		# 需要报证Del/Amp 70%的百分比
		print "#gene-level# $gene info:\n";
		my @type = qw/Amp Del/;
		
		my $gene_info_href = $gene_level_gene_target_idx{$gene}; # ratio/sampleQC/targetQC/CNV_state
		my @target_num = keys %{$gene_info_href};
		
		# skip gene <= 3 targets
		if (@target_num <= 3){
			my $n = @target_num;
			print "gene has $n targets (<=3), will be skipped\n";
			next;
		}

		print(Dumper($gene_info_href));

		for my $t (@type){
			my $ok_idx_aref = &get_cnv_idx($gene_info_href,$gene,$t);
			print "$gene $t (gene-level) ok index is: @{$ok_idx_aref}\n";

			my $seg_aref = &merge_single_exon_into_large_seg($ok_idx_aref);
			print "$gene $t (gene-level) merged single exon seg is: @{$seg_aref}\n";
			
			next if (@{$seg_aref} == 1 and $seg_aref->[0] eq "NA"); # no idx

			# 计算Del/Amp外显子数
			my $cnv_target_num = &get_cnv_target_num($seg_aref);
			my $amp_del_pct = sprintf "%.2f", $cnv_target_num / $exon_num;

			print "$gene has $exon_num targets, $t target num is $cnv_target_num, $t pct is $amp_del_pct\n\n";

			if ($amp_del_pct >= $gene_cutoff){
				# 报出
				my $CN;
				my $ratio_sum;
				my $n;
				my $CNType; # hom-del/het-del/amp
				
				my $ratio_aref = &get_ratio($seg_aref,$gene_info_href);
				
				for my $r (@{$ratio_aref}){
					$n += 1;
					$ratio_sum += $r;
				}

				my $avg_ratio = sprintf "%.2f", $ratio_sum/$n;
				$CN = sprintf "%.1f", $avg_ratio * 2;
				$CNType = &define_CN_state($CN);
				print O2 "$sampleName\t$gene\t$chr\t$start\t$end\t$CN\t$CNType\t$NM\n";
			}else{
				# not gene-level Amp/Del
				next;
			}
		}
	}
}

# RGEM1 promoter dup
#my ($GREM1_dup_state,$copyNumber) = &promoter_dup($infile);
#if ($GREM1_dup_state eq "Amp"){
#	my $cnv_state = &define_CN_state($copyNumber);
#	my $nm = $NM_href->{"GREM1"};
#	print O2 "$sampleName\tGREM1_promoter\tchr15\t32971966\t33004635\t$copyNumber\t$cnv_state\t$nm\n";
#}

close O1;
close O2;


########### SUB FUNCTION DEFINE ###########

sub define_CN_state{
	my ($val) = @_; # CN
	my $state;
	if ($val <= 0.2){
		$state = "hom-del";
	}elsif ($val > 0.2 and $val <= 1.35){
		$state = "het-del";
	}else{
		$state = "amp";
	}

	return($state);
}


sub get_cnv_target_num{
	my ($seg_aref) = @_; # maybe NA
	my $num; # 统计外显子数
	if (@{$seg_aref} == 1 and $seg_aref->[0] eq "NA"){
		# NA
		$num = 0;
	}else{
		for my $seg (@{$seg_aref}){
			# 循环每个大片段
			if ($seg =~ /\_/){
				my @arr = split /\_/, $seg;
				for (@arr){
					$num += 1;
				}
			}else{
				$num += 1;
			}
		}
	}
	return($num);
}

sub get_ratio{
	my ($seg_aref,$gene_href) = @_;
	my @ratio;
	for my $seg (@{$seg_aref}){
		if ($seg =~ /\_/){
			my @arr = split /\_/, $seg;
			for my $i (@arr){
				# for each idx
				my $val = $gene_href->{$i}; # ratio/sampleQC/targetQC/CNV_state
				my $r = (split /\t/, $val)[0];
				next if ($r eq "NA");
				push @ratio, $r;
			}
		}else{
			my $val = $gene_href->{$seg};
			push @ratio, $val;
		}
	}

	return(\@ratio);
}

sub get_cnv_idx{
	my ($gene_info_href,$gene,$cnv_type) = @_;
	my @idx_order = sort {$a <=> $b} (keys %{$gene_info_href});
	
	my @ok_indx; # 返回合格的外显子编号
	for my $idx (@idx_order){
		my $val = $gene_info_href->{$idx};
		my @val = split /\t/, $val; # ratio/sampleQC/targetQC/CNV_state

		my $ratio = $val[0];
		my $sampleQC = $val[1];
		my $targetQC = $val[2];
		my $CNV_state = $val[3];

		next if ($sampleQC eq "Fail");
		next if ($ratio eq "NA" || $targetQC eq "NA" || $targetQC eq "Fail");
		next if ($CNV_state eq "NA");

		if ($cnv_type eq $CNV_state){
			push @ok_indx, $idx;
		}
	}

	if (@ok_indx == 0){
		push @ok_indx, "NA";
	}
	
	return(\@ok_indx);
}


sub merge_single_exon_into_large_seg{
	my ($idx_aref) = @_;

	my @idx = @{$idx_aref};
	my @seg;

	if ($idx[0] eq "NA"){
		# 不存在Del/Amp,返回NA
		push @seg, "NA";
	}elsif (scalar(@idx) == 1 and $idx[0] ne "NA"){
		# 只有1个不用合并
		push @seg, $idx[0];
	}else{
		# >=2个外显子,合并
		# 从1号外显子开始合并. for example: idx is 1/2/3/5/6/8, then @seg will be 1_2_3,5_6,8
		my $first = shift @idx;
		push @seg, $first;

		for my $i (@idx){
			my $last_seg = $seg[-1];
			if ($last_seg =~ /\_/){
				# 之前的seg是合并后的
				my $last_item = (split /\_/, $last_seg)[-1];
				if ($i - $last_item == 1){
					# 相邻合并. first remove last item from @seg, then append current signle value to this ele and use this new seg to replace former seg.
					pop @seg;
					my $append_value = $last_seg."_".$i;
					push @seg, $append_value;
				}else{
					push @seg, $i;
				}
			}else{
				# last item in @seg is a single value
				if ($i - $last_seg == 1){
					# can merge
					pop @seg; # 删除最后一个元素
					my $append_value = $last_seg."_".$i; # 新元素
					push @seg, $append_value; # 替换
				}else{
					# can not merge
					push @seg, $i;
				}
			}
		}
	}
	return(\@seg); # maybe NA
}


sub get_main_NM{
	my ($file) = @_;
	my %main_NM;
	open IN, "$file" or die;
	while (<IN>){
		chomp;
		my @arr = split /\t/;
		$main_NM{$arr[0]} = $arr[1]; # gene=>NM
	}
	close IN;

	# overwrite pre-defined NM
	$main_NM{"BRCA1"} = "NM_007294";
	$main_NM{"BRCA2"} = "NM_000059";
	$main_NM{"MSH2"} = "NM_000251";
	$main_NM{"MLH1"} = "NM_000249";
	$main_NM{"PMS2"} = "NM_000535";
	$main_NM{"MSH6"} = "NM_000179";
	$main_NM{"MUTYH"} = "NM_001128425";
	$main_NM{"EPCAM"} = "NM_002354";

	return(\%main_NM);
}

sub promoter_dup{
	# GREM1基因启动子区是否存在DUP
	my ($file) = @_;
	my @cnv_state;
	my @cn;

	open IN, "$file" or die;
	<IN>;
	while (<IN>){
		chomp;
		my @arr = split /\t/;
		if ($arr[4] eq "SCG5"){
			if ($arr[2] >= 32971966 and $arr[3] <= 32989299){
				push @cnv_state, $arr[-4];
				push @cn, $arr[9]*2;
			}
		}
		if ($arr[4] eq "GREM1"){
			if ($arr[2] >= 32993110 and $arr[3] <= 33004635){
				push @cnv_state, $arr[-4];
				push @cn, $arr[9]*2;
			}
		}
	}
	close IN;

	my $promoter_target_num = scalar @cnv_state; # 9
	my $amp_target_num = 0;
	for my $s (@cnv_state){
		if ($s eq "Amp"){
			$amp_target_num += 1;
		}
	}

	my $GREM1_promoter_cnv_state;
	if ($amp_target_num >= 7){ # 7/9=0.7777778;6/9=0.6666667
		$GREM1_promoter_cnv_state = "Amp";
	}else{
		$GREM1_promoter_cnv_state = "NA";
	}

	my $cn_sum = 0;
	for (@cn){
		$cn_sum += $_;
	}

	my $cn = sprintf "%.2f", $cn_sum/$promoter_target_num;
	print "GREM1_promoter has 9 targets, Amp target num is $amp_target_num\n";
	return($GREM1_promoter_cnv_state,$cn);
}
