use strict;
my $fa = shift;	#'/mnt/isilon/xing_lab/aspera/Feng/20190420feng_293_MIP_exo1_cleanup_T4/extract_unit_splints_consensus_exo1_cleanup_T4.fa'
my $blast = shift; #'/mnt/isilon/xing_lab/aspera/Feng/20190420feng_293_MIP_exo1_cleanup_T4/extract_unit_splints_consensus_exo1_cleanup_T4.blast'
my $splint_length = shift;
my $baseline_umi_f = shift;
my $baseline_umi_r = shift;

my $pre_ID;
my %match_info;
my $min_length_unit = 300;
my $cutoff_lengthfold = .6;
my $ind_umi1 = shift;	#'GGGAA'
my $ind_umi2 = shift;	#'ACTAT'
my %umi_error_range_end = ($ind_umi1 => [1,5], $ind_umi2 => [1,5]);
#my $baseline_umi_f = 20;
#my $baseline_umi_r = 0;
my ($umi1_exist, $umi2_exist) = (0, 0);
if($baseline_umi_f > 0) {
	$umi1_exist = 1;
}
if($baseline_umi_r > 0) {
	$umi2_exist = 1;
}
my $unit_ahead_nt = 10;
my @adjust_umi_order = (0, -1, -2, 1, -3, -4, 2);
my @base_4plus1 = ('A','T','C','G','');
my %umi_possible_seqs_sort;

for my $umi($ind_umi1, $ind_umi2) {
	print "$umi:\n";
	my $umi_possible_seqs_ref = &seq_error1($umi, $umi_error_range_end{$umi}[0], $umi_error_range_end{$umi}[1]);
	my @seq_sort = sort {length($b) <=> length($a) or $a cmp $b} keys %{$umi_possible_seqs_ref};
	$umi_possible_seqs_sort{$umi} = \@seq_sort;
	for my $seq (@seq_sort) {
		print "$seq\t${$umi_possible_seqs_ref}{$seq}\n";
	}
}

my $type;
if (rindex($blast, 'blast') > rindex($blast, '.')) {
	$type = 'blast';
} elsif (rindex($blast, 'hmmer') > rindex($blast, '.')) {
	$type = 'hmmer';
} else {
	die "Input type cannot be identified.";
}

open BLAST, "<", $blast or die "cannot open $blast: $!";
while(<BLAST>) {
	s/\r\n//;
	chomp;
	my @line = split /\s+/;
	next if $line[0] eq '#';
	if ( ($type eq 'blast' and $line[3] >= $splint_length*$cutoff_lengthfold) or ($type eq 'hmmer' and $line[5]-$line[4]+1 >= $splint_length*$cutoff_lengthfold) ) {
		my @info = split '_', $line[0];
		if (exists $match_info{$info[0]}){
			$match_info{$info[0]} = 'NA';
			next;
		} else {
			my $copy_number;
			if (@info == 2) {
				$copy_number = $info[1];
			} else {
				$copy_number = 1;
			}
			if ($type eq 'blast') {
				if ($line[8] > $line[9]) {
					$match_info{$info[0]} = [ '-', $copy_number, $line[2], $line[3], $line[6], $line[7], $line[8], $line[9] ];
				} else {
					$match_info{$info[0]} = [ '+', $copy_number, $line[2], $line[3], $line[6], $line[7], $line[8], $line[9] ];
				}
			} else {
				if($line[6] > $line[7]){
					$match_info{$info[0]} = [ '-', $copy_number, $line[13], $line[5]-$line[4]+1, $line[7], $line[6], $line[5], $line[4] ];
				} else {
					$match_info{$info[0]} = [ '+', $copy_number, $line[13], $line[5]-$line[4]+1, $line[6], $line[7], $line[4], $line[5] ];
				}
			}
		}
	}
}

my $seq;
my $pre_title;
open FA, "<", $fa or die "cannot open $fa: $!";
open OUT, ">", $fa.'.cDNA' or die "cannot write $fa.cDNA: $!";
while(<FA>) {
	s/\r\n//;
	chomp;
	if (/^>/) {
		if (defined $seq and length($seq) >= $min_length_unit){
			my @line = split '_', $pre_title;
			my $readID = substr($line[0], 1, length($line[0])-1);
			if (exists $match_info{$readID} and $match_info{$readID} ne 'NA') {
				my ($strand, $copy_number, $identity, $mapl, $start, $end, $splint_start, $splint_end) = @{$match_info{$readID}};
				#my $seq_UMI_cDNA;
				#my ($best_index1, $best_index2);
				my %adj_pos_seq;
				if ($strand eq '+') {
					my $seq_following = substr( $seq, $end + $splint_length - $splint_end );
					#print "sequence: $seq_following\n" if $readID eq '417601de-1cd8-4dbb-b0b0-b60e84faf39f';
					#$seq_UMI_cDNA = $seq_following;
					for my $i (@adjust_umi_order) {
						my $index_UMI1 = substr($seq_following, $baseline_umi_f+$i, length($ind_umi1));
						my $index_UMI2 = substr($seq_following, length($seq_following)-$unit_ahead_nt-1-$i, length($ind_umi2)); #-length($ind_umi2)  ### remove -1 when no A
						for my $seq1 (@{$umi_possible_seqs_sort{$ind_umi1}}) {
							if (index($index_UMI1, $seq1)>=0) {
								if (!exists $adj_pos_seq{'1'} or ( length($seq1) > length($adj_pos_seq{'1'}[1]) )){
									$adj_pos_seq{'1'} = [$i + index($index_UMI1, $seq1), $seq1];
								}
							}
						}
						for my $seq2 (@{$umi_possible_seqs_sort{$ind_umi2}}) {
							if (index($index_UMI2, $seq2)>=0) {
								if (!exists $adj_pos_seq{'2'} or ( length($seq2) > length($adj_pos_seq{'2'}[1]) )){
									$adj_pos_seq{'2'} = [0-$i + index($index_UMI2, $seq2), $seq2];
								}
							}
						}
					}


					if (exists $adj_pos_seq{'1'} and exists $adj_pos_seq{'2'}) {
						my $pos_umi1_ind = $baseline_umi_f+$adj_pos_seq{'1'}[0];
						my $pos_umi2_ind = length($seq_following)-$baseline_umi_r-$unit_ahead_nt-length($ind_umi2)+$adj_pos_seq{'2'}[0];
						if($umi1_exist == 1) {
							print OUT ">${readID}_UMI1 $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1];
							print OUT "\n", substr($seq_following, 0, $pos_umi1_ind), "\n";
						}
						if($umi2_exist == 1) {
							print OUT ">${readID}_UMI2 $strand $identity $copy_number ", $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
							print OUT "\n", substr( $seq_following, $pos_umi2_ind+length($adj_pos_seq{'2'}[1]), length($seq_following)-($pos_umi2_ind+length($adj_pos_seq{'2'}[1]))-$unit_ahead_nt ), "\n";
						}
						print OUT ">${readID}_cDNA $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1], ' ', $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
						print OUT "\n", substr($seq_following, $pos_umi1_ind+length($adj_pos_seq{'1'}[1]), $pos_umi2_ind-($pos_umi1_ind+length($adj_pos_seq{'1'}[1])) ), "\n";
					}
					

				} else {
					my $seq_following = substr( $seq, $end + $splint_end - $unit_ahead_nt ); ### add -1 when no A
					#$seq_UMI_cDNA = comp_rev($seq_following);

					for my $i (@adjust_umi_order) {
						my $index_UMI1 = substr(comp_rev($seq_following), $baseline_umi_f+$unit_ahead_nt+$i, length($ind_umi1));
						my $index_UMI2 = substr(comp_rev($seq_following), length($seq_following)-$unit_ahead_nt-$i, length($ind_umi2)); #$baseline_umi_r-length($ind_umi2).
						for my $seq1 (@{$umi_possible_seqs_sort{$ind_umi1}}) {
							if (index($index_UMI1, $seq1)>=0) {
								if (!exists $adj_pos_seq{'1'} or ( length($seq1) > length($adj_pos_seq{'1'}[1]) )){
									$adj_pos_seq{'1'} = [$i + index($index_UMI1, $seq1), $seq1];
								}
							} 
						}
						for my $seq2 (@{$umi_possible_seqs_sort{$ind_umi2}}) {
							#print "$index_UMI2\t$seq2:\n" if $readID eq 'c75a3865-2e1b-4838-b430-7c02b7b851e0';
							if (index($index_UMI2, $seq2)>=0) {
								#print index($index_UMI2, $seq2),"\n" if $readID eq 'c75a3865-2e1b-4838-b430-7c02b7b851e0';
								if (!exists $adj_pos_seq{'2'} or ( length($seq2) > length($adj_pos_seq{'2'}[1]) )){
									$adj_pos_seq{'2'} = [0-$i + index($index_UMI2, $seq2), $seq2];
								}
							}
						}
					}

					if (exists $adj_pos_seq{'1'} and exists $adj_pos_seq{'2'}) {
						my $pos_umi1_ind = $baseline_umi_f+$unit_ahead_nt+$adj_pos_seq{'1'}[0];
						my $pos_umi2_ind = length($seq_following)-$baseline_umi_r-length($ind_umi2)+$adj_pos_seq{'2'}[0];
						if($umi1_exist == 1) {
							print OUT ">${readID}_UMI1 $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1];
							print OUT "\n", substr(comp_rev($seq_following), $unit_ahead_nt, $pos_umi1_ind-$unit_ahead_nt), "\n";
						}
						if($umi2_exist == 1) {
							print OUT ">${readID}_UMI2 $strand $identity $copy_number ", $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
							print OUT "\n", substr(comp_rev($seq_following), $pos_umi2_ind+length($adj_pos_seq{'2'}[1])), "\n";
						}
						print OUT ">${readID}_cDNA $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1], ' ', $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
						print OUT "\n", substr(comp_rev($seq_following), $pos_umi1_ind+length($adj_pos_seq{'1'}[1]), $pos_umi2_ind-($pos_umi1_ind+length($adj_pos_seq{'1'}[1])) ), "\n";
					}
	
				}
				
				print "$readID";
				print "\t$_" for @{$match_info{$readID}};
				print "\t", length($seq);

				if (exists $adj_pos_seq{'1'}) {
					print "\t",$adj_pos_seq{'1'}[0],"\t",$adj_pos_seq{'1'}[1];
				} else {
					print "\tNA\tNA";
				}
				if (exists $adj_pos_seq{'2'}) {
					print "\t",$adj_pos_seq{'2'}[0],"\t",$adj_pos_seq{'2'}[1];
				} else {
					print "\tNA\tNA";
				}
				print "\n";
			}
		}
		$pre_title = $_;
		$seq = '';
	} else {
		$seq .= $_;
	}
}

if (defined $seq and length($seq) >= $min_length_unit){
			my @line = split '_', $pre_title;
			my $readID = substr($line[0], 1, length($line[0])-1);
			if (exists $match_info{$readID} and $match_info{$readID} ne 'NA') {
				my ($strand, $copy_number, $identity, $mapl, $start, $end, $splint_start, $splint_end) = @{$match_info{$readID}};
				#my $seq_UMI_cDNA;
				#my ($best_index1, $best_index2);
				my %adj_pos_seq;
				if ($strand eq '+') {
					my $seq_following = substr( $seq, $end + $splint_length - $splint_end );
					#print "sequence: $seq_following\n" if $readID eq '417601de-1cd8-4dbb-b0b0-b60e84faf39f';
					#$seq_UMI_cDNA = $seq_following;
					for my $i (@adjust_umi_order) {
						my $index_UMI1 = substr($seq_following, $baseline_umi_f+$i, length($ind_umi1));
						my $index_UMI2 = substr($seq_following, length($seq_following)-$unit_ahead_nt-1-$i, length($ind_umi2)); #-length($ind_umi2)  ### remove -1 when no A
						for my $seq1 (@{$umi_possible_seqs_sort{$ind_umi1}}) {
							if (index($index_UMI1, $seq1)>=0) {
								if (!exists $adj_pos_seq{'1'} or ( length($seq1) > length($adj_pos_seq{'1'}[1]) )){
									$adj_pos_seq{'1'} = [$i + index($index_UMI1, $seq1), $seq1];
								}
							}
						}
						for my $seq2 (@{$umi_possible_seqs_sort{$ind_umi2}}) {
							if (index($index_UMI2, $seq2)>=0) {
								if (!exists $adj_pos_seq{'2'} or ( length($seq2) > length($adj_pos_seq{'2'}[1]) )){
									$adj_pos_seq{'2'} = [0-$i + index($index_UMI2, $seq2), $seq2];
								}
							}
						}
					}


					if (exists $adj_pos_seq{'1'} and exists $adj_pos_seq{'2'}) {
						my $pos_umi1_ind = $baseline_umi_f+$adj_pos_seq{'1'}[0];
						my $pos_umi2_ind = length($seq_following)-$baseline_umi_r-$unit_ahead_nt-length($ind_umi2)+$adj_pos_seq{'2'}[0];
						if($umi1_exist == 1) {
							print OUT ">${readID}_UMI1 $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1];
							print OUT "\n", substr($seq_following, 0, $pos_umi1_ind), "\n";
						}
						if($umi2_exist == 1) {
							print OUT ">${readID}_UMI2 $strand $identity $copy_number ", $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
							print OUT "\n", substr( $seq_following, $pos_umi2_ind+length($adj_pos_seq{'2'}[1]), length($seq_following)-($pos_umi2_ind+length($adj_pos_seq{'2'}[1]))-$unit_ahead_nt ), "\n";
						}
						print OUT ">${readID}_cDNA $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1], ' ', $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
						print OUT "\n", substr($seq_following, $pos_umi1_ind+length($adj_pos_seq{'1'}[1]), $pos_umi2_ind-($pos_umi1_ind+length($adj_pos_seq{'1'}[1])) ), "\n";
					}
					

				} else {
					my $seq_following = substr( $seq, $end + $splint_end - $unit_ahead_nt ); ### add -1 when no A
					#$seq_UMI_cDNA = comp_rev($seq_following);

					for my $i (@adjust_umi_order) {
						my $index_UMI1 = substr(comp_rev($seq_following), $baseline_umi_f+$unit_ahead_nt+$i, length($ind_umi1));
						my $index_UMI2 = substr(comp_rev($seq_following), length($seq_following)-$unit_ahead_nt-$i, length($ind_umi2)); #$baseline_umi_r-length($ind_umi2).
						for my $seq1 (@{$umi_possible_seqs_sort{$ind_umi1}}) {
							if (index($index_UMI1, $seq1)>=0) {
								if (!exists $adj_pos_seq{'1'} or ( length($seq1) > length($adj_pos_seq{'1'}[1]) )){
									$adj_pos_seq{'1'} = [$i + index($index_UMI1, $seq1), $seq1];
								}
							} 
						}
						for my $seq2 (@{$umi_possible_seqs_sort{$ind_umi2}}) {
							#print "$index_UMI2\t$seq2:\n" if $readID eq 'c75a3865-2e1b-4838-b430-7c02b7b851e0';
							if (index($index_UMI2, $seq2)>=0) {
								#print index($index_UMI2, $seq2),"\n" if $readID eq 'c75a3865-2e1b-4838-b430-7c02b7b851e0';
								if (!exists $adj_pos_seq{'2'} or ( length($seq2) > length($adj_pos_seq{'2'}[1]) )){
									$adj_pos_seq{'2'} = [0-$i + index($index_UMI2, $seq2), $seq2];
								}
							}
						}
					}

					if (exists $adj_pos_seq{'1'} and exists $adj_pos_seq{'2'}) {
						my $pos_umi1_ind = $baseline_umi_f+$unit_ahead_nt+$adj_pos_seq{'1'}[0];
						my $pos_umi2_ind = length($seq_following)-$baseline_umi_r-length($ind_umi2)+$adj_pos_seq{'2'}[0];
						if($umi1_exist == 1) {
							print OUT ">${readID}_UMI1 $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1];
							print OUT "\n", substr(comp_rev($seq_following), $unit_ahead_nt, $pos_umi1_ind-$unit_ahead_nt), "\n";
						}
						if($umi2_exist == 1) {
							print OUT ">${readID}_UMI2 $strand $identity $copy_number ", $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
							print OUT "\n", substr(comp_rev($seq_following), $pos_umi2_ind+length($adj_pos_seq{'2'}[1])), "\n";
						}
						print OUT ">${readID}_cDNA $strand $identity $copy_number ", $adj_pos_seq{'1'}[0].':'.$adj_pos_seq{'1'}[1], ' ', $adj_pos_seq{'2'}[0].':'.$adj_pos_seq{'2'}[1];
						print OUT "\n", substr(comp_rev($seq_following), $pos_umi1_ind+length($adj_pos_seq{'1'}[1]), $pos_umi2_ind-($pos_umi1_ind+length($adj_pos_seq{'1'}[1])) ), "\n";
					}
	
				}
				
				print "$readID";
				print "\t$_" for @{$match_info{$readID}};
				print "\t", length($seq);

				if (exists $adj_pos_seq{'1'}) {
					print "\t",$adj_pos_seq{'1'}[0],"\t",$adj_pos_seq{'1'}[1];
				} else {
					print "\tNA\tNA";
				}
				if (exists $adj_pos_seq{'2'}) {
					print "\t",$adj_pos_seq{'2'}[0],"\t",$adj_pos_seq{'2'}[1];
				} else {
					print "\tNA\tNA";
				}
				print "\n";
			}
}


	sub index_all {
		my ($seq, $seed) = @_;
		my @index;
		my $initial = -1;
		while(1) {
			$initial = index($seq, $seed, $initial+1);
			#print "$initial\n";
			if ($initial >= 0) {
				push @index, $initial;
			} else {
				last;
			}
		}
		\@index;
	}


	sub comp_rev {
		my $seq = reverse($_[0]);
		$seq =~ tr/ATCG/TAGC/;
		$seq;
	}




sub seq_error1 {
	my ($seq, $range_start, $range_end) = @_;
	my @bases = split '', $seq;
	my %error_seqs;
	for my $i ($range_start .. $range_end) {
		my @bases_tmp = @bases;
		for my $n (@base_4plus1) {
			$bases_tmp[$i-1] = $n;
			my $seq_tmp = join '', @bases_tmp;
			$error_seqs{$seq_tmp} ++;
		}
	}
	\%error_seqs;
}




