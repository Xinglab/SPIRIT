use strict;
my $fa = shift;	#'/mnt/isilon/xing_lab/aspera/Feng/20190503feng_RT_TS_noUMI2/RT_TS_noUMI2.fa'
my $blast = shift; #'/mnt/isilon/xing_lab/aspera/Feng/20190503feng_RT_TS_noUMI2/RT_TS_noUMI2_T8_max.hmmer1_sort'
my $splint_length = shift; #'393'
my $min_matched_splint_length = 30;
my $cutoff_splint_lengthfold = 1.5;
my $cutoff_unit_length_ratio = 0.85;

my $pre_ID;
my %align_start_end;
my %align_lengths;
my %seqs_between_splint;

my $type;
if (rindex($blast, 'blast') > rindex($blast, '.')) {
	$type = 'blast';
} elsif (rindex($blast, 'hmmer') > rindex($blast, '.')) {
	$type = 'hmmer';
} else {
	die "Input type cannot be identified.";
}

open BLAST, "<", $blast or die "cannot open $blast: $!";
open READ_INFO, ">", $blast.'.info' or die "cannot write $blast.info: $!";
while(<BLAST>) {
	s/\r\n//;
	chomp;
	my @line = split /\s+/;
	next if $line[0] eq '#';
	if ( $line[0] ne $pre_ID and $align_lengths{$pre_ID} >= $splint_length*$cutoff_splint_lengthfold and exists $align_start_end{$pre_ID}) {
		my @sort_start_end = sort {${$a}[0] <=> ${$b}[0] or ${$a}[1] <=> ${$b}[1]} @{$align_start_end{$pre_ID}};
		my $between_info_ref = &cDNAs_in_raw_windows(@sort_start_end);
		my @between_info = @{$between_info_ref};

		printf READ_INFO "$pre_ID\t%.3f", $align_lengths{$pre_ID}/$splint_length;
		print READ_INFO "\t$between_info[$_]" for (0 .. 3);

		if ($between_info[-1] ne 'NA'){
			$seqs_between_splint{$pre_ID} = $between_info_ref;
			print READ_INFO "\t", scalar(@{$between_info[-1]}),"\t", $between_info[-1][0][2], "\t", $between_info[-1][-1][2], "\n";
			print READ_INFO "unit_length:";
			print READ_INFO "\t${$_}[2]" for @{$between_info[-1]};
			print READ_INFO "\n";
			print READ_INFO "mean_score:";
			print READ_INFO "\t${$_}[3]" for @{$between_info[-1]};
			print READ_INFO "\n";
		} else {
			print READ_INFO "\t0\n";
		}

	} elsif (defined($pre_ID) and $line[0] ne $pre_ID) {
		printf READ_INFO "$pre_ID\t%.3f\tNA\n", $align_lengths{$pre_ID}/$splint_length;
	}
	if ($type eq 'blast') {
		push @{$align_start_end{$line[0]}}, [$line[6], $line[7], $line[8], $line[9], $line[11]] if $line[3]>=$min_matched_splint_length;
		$align_lengths{$line[0]} += $line[3];
	} else {
		if($line[6]<$line[7]){
			push @{$align_start_end{$line[0]}}, [$line[6], $line[7], $line[4], $line[5], $line[13]] if $line[5]-$line[4]+1>=$min_matched_splint_length;
		} else {
			push @{$align_start_end{$line[0]}}, [$line[7], $line[6], $line[5], $line[4], $line[13]] if $line[5]-$line[4]+1>=$min_matched_splint_length;
		}
		$align_lengths{$line[0]} += $line[5]-$line[4]+1;
	}
	
	$pre_ID = $line[0];
}
if ( $align_lengths{$pre_ID} >= $splint_length*$cutoff_splint_lengthfold and exists $align_start_end{$pre_ID}) {
		my @sort_start_end = sort {${$a}[0] <=> ${$b}[0] or ${$a}[1] <=> ${$b}[1]} @{$align_start_end{$pre_ID}};
		my $between_info_ref = &cDNAs_in_raw_windows(@sort_start_end);
		my @between_info = @{$between_info_ref};

		printf READ_INFO "$pre_ID\t%.3f", $align_lengths{$pre_ID}/$splint_length;
		print READ_INFO "\t$between_info[$_]" for (0 .. 3);

		if ($between_info[-1] ne 'NA'){
			$seqs_between_splint{$pre_ID} = $between_info_ref;
			print READ_INFO "\t", scalar(@{$between_info[-1]}),"\t", $between_info[-1][0][2], "\t", $between_info[-1][-1][2], "\n";
			print READ_INFO "unit_length:";
			print READ_INFO "\t${$_}[2]" for @{$between_info[-1]};
			print READ_INFO "\n";
			print READ_INFO "mean_score:";
			print READ_INFO "\t${$_}[3]" for @{$between_info[-1]};
			print READ_INFO "\n";
		} else {
			print READ_INFO "\t0\n";
		}

} else {
	printf READ_INFO "$pre_ID\t%.3f\tNA\n", $align_lengths{$pre_ID}/$splint_length;
}

my $seq;
my $pre_title;
open FA, "<", $fa or die "cannot open $fa: $!";
while(<FA>) {
	s/\r\n//;
	chomp;
	if (/^>/) {
		if (defined $seq){
			my @line = split /\s+/, $pre_title;
			my $readID = substr($line[0], 1, length($line[0])-1);
			if (exists $seqs_between_splint{$readID}) {
				my ($raw_total, $major_strand, $between_2_ends, $strand, $info_ref) = @{$seqs_between_splint{$readID}};
				for my $seq_coord (@{$info_ref}) {

					my ($start, $end, $length, $mean_identity) = @{$seq_coord};
					if ($start + $length > length($seq)){
						$end = length($seq);
						$length = length($seq) - $start;
					} if ($start < 0){
						$length = $length + $start;
						$start = 0;
					}
					my $between_splints_seq = substr($seq, $start, $length);
					print ">${readID}_${start}_${end}_${length}_${mean_identity}\n$between_splints_seq\n"; #_${raw_total}_${major_strand}_${between_2_ends}_${strand}
				}
			}
		}
		$pre_title = $_;
		$seq = '';
	} else {
		$seq .= $_;
	}
}

if (defined $seq){
			my @line = split /\s+/, $pre_title;
			my $readID = substr($line[0], 1, length($line[0])-1);
			if (exists $seqs_between_splint{$readID}) {
				my ($raw_total, $major_strand, $between_2_ends, $strand, $info_ref) = @{$seqs_between_splint{$readID}};
				for my $seq_coord (@{$info_ref}) {

					my ($start, $end, $length, $mean_identity) = @{$seq_coord};
					if ($start + $length > length($seq)){
						$end = length($seq);
						$length = length($seq) - $start;
					} if ($start < 0){
						$length = $length + $start;
						$start = 0;
					}
					my $between_splints_seq = substr($seq, $start, $length);
					print ">${readID}_${start}_${end}_${length}_${mean_identity}\n$between_splints_seq\n"; #_${raw_total}_${major_strand}_${between_2_ends}_${strand}
				}
			}
}

#print scalar(keys %seqs_between_splint), "\n";

sub cDNAs_in_raw_windows{
	my @sorted_2ends = @_;
	my $is_different_strand = 1;
	my @mt0_pairs = grep {${$_}[3]-${$_}[2]>0} @sorted_2ends;
	$is_different_strand = 0 if @mt0_pairs == @sorted_2ends or @mt0_pairs == 0;
	my @major_strand_pairs = ();
	my @adj_pairs;
	my $strand = 'NA';
	if (@mt0_pairs >= @sorted_2ends/2 and @mt0_pairs >= 2) {
		$strand = '+';
		@major_strand_pairs = @mt0_pairs;
		@adj_pairs = map { [ ${$_}[2]-1, $splint_length-${$_}[3] ] } @mt0_pairs;
	} elsif (@mt0_pairs <= @sorted_2ends/2 and @sorted_2ends-@mt0_pairs >= 2) {
		$strand = '-';
		@major_strand_pairs = grep {${$_}[3]-${$_}[2]<0} @sorted_2ends;
		@adj_pairs = map { [ $splint_length-${$_}[2], ${$_}[3]-1 ] } @major_strand_pairs;
	}
	my @sorted_2ends_adjusted;
	for my $i (0 .. $#major_strand_pairs) {
		push @sorted_2ends_adjusted, [ $major_strand_pairs[$i][0]-$adj_pairs[$i][0], $major_strand_pairs[$i][1]+$adj_pairs[$i][1], $major_strand_pairs[$i][4]];
	}
	my @between_2ends_adjusted = ();
	my $total_length = 0;
	for my $i (0 .. $#sorted_2ends_adjusted-1) {
		if ($sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][1]>0){
			push @between_2ends_adjusted, [ $sorted_2ends_adjusted[$i][0] - 10, $sorted_2ends_adjusted[$i+1][0]-1 + 10, $sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][0]+20, ($sorted_2ends_adjusted[$i+1][2]+$sorted_2ends_adjusted[$i][2])*.5 ];
			$total_length += $sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][0]+20;
		}
		
	}
	if (@between_2ends_adjusted == 0) {
		return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, 'NA'];
	} elsif (@between_2ends_adjusted == 1) {
		return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, \@between_2ends_adjusted];
	} elsif (@between_2ends_adjusted == 2) {
		my @sort_between_2ends_adjusted = sort {${$a}[2] <=> ${$b}[2]} @between_2ends_adjusted;
		if ($sort_between_2ends_adjusted[0][2] > $sort_between_2ends_adjusted[-1][2]*$cutoff_unit_length_ratio) {
			my @sort_between_2ends_adjusted2 = sort {${$a}[3] <=> ${$b}[3]} @sort_between_2ends_adjusted;
			return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, [$sort_between_2ends_adjusted2[0]]];
		} else {
			return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, [$sort_between_2ends_adjusted[0]]];
		}
	} else {
		my @sort_between_2ends_adjusted = sort {${$a}[2] <=> ${$b}[2]} @between_2ends_adjusted;
		my @sort_between_2ends_adjusted_group;
		for my $i (0 .. $#sort_between_2ends_adjusted) {
			push @{$sort_between_2ends_adjusted_group[$i]}, $sort_between_2ends_adjusted[$i];
			for my $j ($i+1 .. $#sort_between_2ends_adjusted) {
				if ($sort_between_2ends_adjusted[$i][2]/$sort_between_2ends_adjusted[$j][2] >= $cutoff_unit_length_ratio) {
					push @{$sort_between_2ends_adjusted_group[$i]}, $sort_between_2ends_adjusted[$j];
				} else {
					last;
				}
			}
		}
		my @sort_between_2ends_adjusted2 = sort { scalar(@{$b}) <=> scalar(@{$a}) or (${$a}[-1][2]-${$a}[0][2]) <=> (${$b}[-1][2]-${$b}[0][2]) } @sort_between_2ends_adjusted_group;
		if (@{$sort_between_2ends_adjusted2[0]} == 2) {
			my @sort_between_2ends_adjusted3 = sort {${$b}[3] <=> ${$a}[3]} @{$sort_between_2ends_adjusted2[0]};
			return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, [$sort_between_2ends_adjusted3[0]]];
		} else {
			return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, $sort_between_2ends_adjusted2[0]];
		}
	}
}

sub cDNAs_in_raw_improved{
	my @sorted_2ends = @_;
	my $is_different_strand = 1;
	my @mt0_pairs = grep {${$_}[3]-${$_}[2]>0} @sorted_2ends;
	$is_different_strand = 0 if @mt0_pairs == @sorted_2ends or @mt0_pairs == 0;
	my @major_strand_pairs = ();
	my @adj_pairs;
	my $strand = 'NA';
	if (@mt0_pairs >= @sorted_2ends/2 and @mt0_pairs >= 2) {
		$strand = '+';
		@major_strand_pairs = @mt0_pairs;
		@adj_pairs = map { [ ${$_}[2]-1, $splint_length-${$_}[3] ] } @mt0_pairs;
	} elsif (@mt0_pairs <= @sorted_2ends/2 and @sorted_2ends-@mt0_pairs >= 2) {
		$strand = '-';
		@major_strand_pairs = grep {${$_}[3]-${$_}[2]<0} @sorted_2ends;
		@adj_pairs = map { [ $splint_length-${$_}[2], ${$_}[3]-1 ] } @major_strand_pairs;
	}
	my @sorted_2ends_adjusted;
	for my $i (0 .. $#major_strand_pairs) {
		push @sorted_2ends_adjusted, [ $major_strand_pairs[$i][0]-$adj_pairs[$i][0], $major_strand_pairs[$i][1]+$adj_pairs[$i][1], $major_strand_pairs[$i][4]];
	}
	my @between_2ends_adjusted = ();
	my $total_length = 0;
	for my $i (0 .. $#sorted_2ends_adjusted-1) {
		if ($sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][1]>0){
			push @between_2ends_adjusted, [ $sorted_2ends_adjusted[$i][0] - 10, $sorted_2ends_adjusted[$i+1][0]-1 + 10, $sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][0]+20, ($sorted_2ends_adjusted[$i+1][2]+$sorted_2ends_adjusted[$i][2])*.5 ];
			$total_length += $sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][0]+20;
		}
		
	}
	if (@between_2ends_adjusted == 0) {
		return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, 'NA'];
	} elsif (@between_2ends_adjusted == 1) {
		return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, \@between_2ends_adjusted];
	} elsif (@between_2ends_adjusted == 2) {
		my @sort_between_2ends_adjusted = sort {${$a}[2] <=> ${$b}[2]} @between_2ends_adjusted;
		if ($sort_between_2ends_adjusted[0][2] > $sort_between_2ends_adjusted[-1][2]*$cutoff_unit_length_ratio) {
			my @sort_between_2ends_adjusted2 = sort {${$b}[3] <=> ${$a}[3]} @sort_between_2ends_adjusted;
			return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, [$sort_between_2ends_adjusted2[0]]];
		} else {
			return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, [$sort_between_2ends_adjusted[0]]];
		}
	} else {
		my @sort_between_2ends_adjusted = sort {${$a}[2] <=> ${$b}[2]} @between_2ends_adjusted;
		while (1) {
			if ($sort_between_2ends_adjusted[0][2] > $sort_between_2ends_adjusted[-1][2]*$cutoff_unit_length_ratio) {
				if (@sort_between_2ends_adjusted == 2) {
					my @sort_between_2ends_adjusted2 = sort {${$b}[3] <=> ${$a}[3]} @sort_between_2ends_adjusted;
					return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, [$sort_between_2ends_adjusted2[0]]];
				} else {
					return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, \@sort_between_2ends_adjusted];
				}
			} elsif (@sort_between_2ends_adjusted == 2) {
				return [scalar(@sorted_2ends), scalar(@major_strand_pairs), scalar(@between_2ends_adjusted), $strand, [$sort_between_2ends_adjusted[0]]];
			} else {
				my $median_length_less1 = $sort_between_2ends_adjusted[int((@sort_between_2ends_adjusted+1)/2-1)][2];
				my $median_length_less2 = $sort_between_2ends_adjusted[int((@sort_between_2ends_adjusted)/2)][2];
				if ($median_length_less1/$sort_between_2ends_adjusted[0][2] <= $sort_between_2ends_adjusted[-1][2]/$median_length_less2) {
					pop @sort_between_2ends_adjusted;
				} else {
					shift @sort_between_2ends_adjusted;
				}
			}
		}
	}
}

sub cDNAs_in_raw_quick{
	my @sorted_2ends = @_;
	my @mt0_pairs = grep {${$_}[3]-${$_}[2]>0} @sorted_2ends;
	my @adj_pairs;
	my $strand;
	if (@mt0_pairs == @sorted_2ends) {
		$strand = '+';
		@adj_pairs = map { [ ${$_}[2]-1, $splint_length-${$_}[3] ] } @sorted_2ends;
	} elsif (@mt0_pairs == 0) {
		$strand = '-';
		@adj_pairs = map { [ $splint_length-${$_}[2], ${$_}[3]-1 ] } @sorted_2ends;
	} else {
		return [0, 'NA'];
	}
	my @sorted_2ends_adjusted;
	for my $i (0 .. $#sorted_2ends) {
		push @sorted_2ends_adjusted, [ $sorted_2ends[$i][0]-$adj_pairs[$i][0], $sorted_2ends[$i][1]+$adj_pairs[$i][1] ];
	}
	my @between_2ends_adjusted;
	my $total_length = 0;
	for my $i (0 .. $#sorted_2ends_adjusted-1) {
		push @between_2ends_adjusted, [ $sorted_2ends_adjusted[$i][0] - 10, $sorted_2ends_adjusted[$i+1][0]-1 + 10, $sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][0]+20 ] if $sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][1]>0;
		$total_length += $sorted_2ends_adjusted[$i+1][0]-1-$sorted_2ends_adjusted[$i][0]+20;
	}
	if ($total_length <= 0 or @between_2ends_adjusted == 0) {
		return [-1, $strand];
	}
	my @sort_between_2ends_adjusted = sort {${$a}[2] <=> ${$b}[2]} @between_2ends_adjusted;
	my $median_length_less = $sort_between_2ends_adjusted[int((@sort_between_2ends_adjusted+1)/2-1)][2];
	#my $isAboutAve = grep { ${$_}[2]/$total_length*@between_2ends_adjusted > .8 and ${$_}[2]/$total_length*@between_2ends_adjusted < 1.25 } @between_2ends_adjusted;
	my @isAboutMedian = grep { ${$_}[2]/$median_length_less > .8 and ${$_}[2]/$median_length_less < 1.25 } @between_2ends_adjusted;
	#if ($isAboutAve == @sorted_2ends_adjusted-1) {
	#	return [1, $strand, \@between_2ends_adjusted];
	#} else {
	#	return [-1, $strand];
	#}
	if(@isAboutMedian >= 1) {
		return [1, $strand, \@isAboutMedian];
	} else {
		return [-1, $strand];
	}
}

