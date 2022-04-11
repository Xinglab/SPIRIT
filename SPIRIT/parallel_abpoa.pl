use strict;
use threads;

my $input = shift;
my ($input_dir, $fa);
if (rindex($input, '/')>=0){
	$input_dir = substr($input, 0, rindex($input, '/')+1);
	#$fa = substr($input, rindex($input, '/')+1);
} else {
	$input_dir = './';
}
$fa = $input;

#$max_thread = 2 if !defined $max_thread or $max_thread<2;
my $out = shift;
my $max_thread = shift;
if(!defined($max_thread)) {
	$max_thread = 1;
}
my %split_files_key_read;
my $quiet = 'no';
my $max_units_length = 100_000;
&split_fa_file($fa, $max_thread);


if ($max_thread >= 2) {
	my @ths;
	while (my ($file, $key_read) = each %split_files_key_read) {
		my $th = threads -> new({'context' => 'void'}, \&parallel_abPOA, [$file, $key_read]);
		my $th_id = $th->tid();
		print " Worker $th_id begins to process $file.\n" if $quiet eq 'no';
		push @ths, $th;
	}
	for (@ths) {
		my $th_id = $_->tid();
		$_ -> join();
		print " Worker $th_id finished reporting.\n" if $quiet eq 'no';
	}
	@ths = ();
} else {
	while (my ($file, $key_read) = each %split_files_key_read) {
		system("rm -f $input_dir$file.pir") if -f "$input_dir$file.pir";
		system("rm -f $input_dir$file.log2") if -f "$input_dir$file.log2";
		&parallel_abPOA([$file, $key_read]);
	}
}

system("rm $out") if -f $out;
while (my ($file, $key_read) = each %split_files_key_read) {
	system("cat $input_dir${file}_single.fa >> $out");
	system("rm $input_dir${file}_single.fa");
	system("cat $input_dir${file}.pir2 >> $out");
	system("rm $input_dir${file}.pir2") if $max_thread >= 2;
	system("rm $input_dir${file}.pir") if $max_thread >= 2;
	system("rm $input_dir${file}") if $max_thread >= 2;
}

	sub parallel_abPOA {
		my ($file, $key_read) = @{$_[0]};
		my ($pre_ID, @reads);
		my $tag = 0;
		my (@failed_files, @succeed_files);
		open IN, "<", $input_dir.$file or die "cannot open $input_dir$file: $!";
		open SINGLEFA, ">", $input_dir.$file.'_single.fa' or die "cannot write $input_dir${file}_single.fa: $!";
		while(<IN>) {
			s/\r\n//;
			chomp;
			if (substr($_, 0, 1) eq '>'){
				my @line = split '_';
				if ( defined $pre_ID and $line[0] ne $pre_ID ) {
					if ($tag == 1) {
						if (@reads/2 >= 2) {
							my $read_fa = join " ", @reads;
							my $command;
							#$read_fa =~ tr/>/~/;
							if (length($read_fa) <= $max_units_length) {
								$command = "echo \'$read_fa\' | tr \" \" \"\n\" | abpoa -r 0 - >> $input_dir$file.pir 2>>$input_dir$file.log2";
							} else {
								my $head_read_fa = substr($read_fa, 0, $max_units_length);
								my $space_rindex = rindex($head_read_fa, " >");
								my $space_index = index($head_read_fa, " >");
								if ($space_rindex > $space_index) {
									$head_read_fa = substr($head_read_fa, 0, $space_rindex);
									$command = "echo \'$head_read_fa\' | tr \" \" \"\n\" | abpoa -r 0 - >> $input_dir$file.pir 2>>$input_dir$file.log2";
								} else {
									print SINGLEFA "$reads[0]\n$reads[1]\n";
								}
							}
							if (defined $command) {
								my $status = system($command);
								my $exit_code = ($status >> 8) & 0xff;
								if ($exit_code != 0) {
									push @failed_files, $pre_ID.':'.$exit_code.':'.$status;
								} else {
									push @succeed_files, $pre_ID.'_'.(@reads/2);
								}
							}
						} else {
							print SINGLEFA "$_\n" for @reads;
						}
					}
					@reads = ();
				}
				$pre_ID = $line[0];
				if ($line[0] eq $key_read){
					$tag = 1;
				}
			}
			push @reads, $_;
		}

		if (@reads/2 >= 2) {
			my $read_fa = join " ", @reads;
			#$read_fa =~ tr/>/~/;
			my $command;
			#$read_fa =~ tr/>/~/;
			if (length($read_fa) <= $max_units_length) {
				$command = "echo \'$read_fa\' | tr \" \" \"\n\" | abpoa -r 0 - >> $input_dir$file.pir 2>>$input_dir$file.log2";
			} else {
				my $head_read_fa = substr($read_fa, 0, $max_units_length);
				my $space_rindex = rindex($head_read_fa, " >");
				my $space_index = index($head_read_fa, " >");
				if ($space_rindex > $space_index) {
					$head_read_fa = substr($head_read_fa, 0, $space_rindex);
					$command = "echo \'$head_read_fa\' | tr \" \" \"\n\" | abpoa -r 0 - >> $input_dir$file.pir 2>>$input_dir$file.log2";
				} else {
					print SINGLEFA "$reads[0]\n$reads[1]\n";
				}
			}
			if (defined $command) {
				my $status = system($command);
				my $exit_code = ($status >> 8) & 0xff;
				if ($exit_code != 0) {
					push @failed_files, $pre_ID.':'.$exit_code.':'.$status;
				} else {
					push @succeed_files, $pre_ID.'_'.(@reads/2);
				}
			}
		} else {
			print SINGLEFA "$_\n" for @reads;
		}

		close IN;
		close SINGLEFA;

		if (-f $input_dir.$file.'.pir') {
			my $n = 0;
			open PIR, "<", $input_dir.$file.'.pir' or die "cannot open $input_dir$file.pir: $!";
			open PIROUT, ">", $input_dir.$file.'.pir2' or die "cannot write $input_dir$file.pir2: $!";
			while (<PIR>) {
				chomp;
				if (substr($_, 0, 1) eq '>') {
					print PIROUT "$succeed_files[$n]\n";
					$n ++;
				} else {
					print PIROUT "$_\n";
				}
			}
			close PIR;
			close PIROUT;
			my $m = @succeed_files;
			die "$n vs. $m! Read IDs may be wrong." if $n != $m;
		}

		if (@failed_files > 0) {
			#record to log file if failing to estimate
			open FAIL, ">>", $input_dir.$file.".log";
			print FAIL "abpoa failed to process the following FASTA files:\n@failed_files\n";
			close FAIL;
		}
	}

	sub split_fa_file {
		my ($in, $n) = @_;
		my $in_raw;
		if (rindex($in, "/")>=0) {
			$in_raw = (substr($in, rindex($in, "/")+1));
		} else {
			$in_raw = $in;
		}
		if ($n >= 2) {
			#$input_dir = $output_dir;
			my $size = -s $in;
			if($size < 1_001){
				print "Fail to get file size for $in.\nFatal error. Aborted.\n";
				die "Fail to get file size for $in: $!";
			}
			my $division;
			if ($size%$n != 0) {
				$division = int($size/$n)+1;
			} else {
				$division = int($size/$n);
			}
			if(!defined($division) or $division < 0.01*$size){
				print "Fail to calculate divided file size, total file size $size, thread $n.\nFatal error. Aborted.\n";
				die "Fail to calculate divided file size, total file size $size, thread $n: $!";
			}
			system "split -b $division $in $input_dir$in_raw";
			my @split_sam = <$input_dir$in_raw*>;
			print " Divided SAM sizes:\n" if $quiet eq 'no';;
			my $n2 = 0;
			for my $sam (@split_sam) {
				if ($sam=~/$in_raw[a-z]+$/){
					$n2 ++;
					printf " $sam\t%15d\n", (-s $sam) if $quiet eq 'no';;
				}
			}
			if ($n2 == $n){
				print " FA was divided successfully.\n" if $quiet eq 'no';;
			} else {
				print "Cannot split $in into $n ($n2) pieces with size of $division and named them as $input_dir$in_raw.\nFatal error. Aborted.\n";
				die "Cannot split $in into $n ($n2) pieces with size of $division and named them as $input_dir$in_raw: $!";
			}
			my @split_files;
			opendir DIR, $input_dir or die "cannot open directory $input_dir: $!";
			for my $file (readdir DIR) {
				if ($file=~/$in_raw[a-z]+$/) {
					push @split_files, $file;
				}
			}
			closedir DIR;
			print " First read of divided FA files: \n" if $quiet eq 'no';;
			my @split_files_sort = sort{$a cmp $b} @split_files;
			for my $i (1 .. $#split_files_sort) {
				my (%read_name, @add_reads);
				open FILE, "<", $input_dir.$split_files_sort[$i] or die "cannot open $input_dir$split_files_sort[$i]: $!";
				while (<FILE>) {
					chomp;
					#my @line = split /\t/;
					if (substr($_, 0, 1) eq '>'){
						my @line = split '_';
						$read_name{$line[0]} ++;
						if (scalar(keys %read_name) >= 2) {
							$split_files_key_read{$split_files_sort[$i]} = $line[0];
							print " $split_files_sort[$i]: $line[0]\n" if $quiet eq 'no';;
							open OUTPUT, ">>", $input_dir.$split_files_sort[$i-1] or die;
							print OUTPUT "$_\n" for @add_reads;
							close OUTPUT;
							@add_reads = ();
							last;
						}
					}
					push @add_reads, $_;
				}
			}
			open FILE, "<", $input_dir.$split_files_sort[0] or die "cannot open $input_dir$split_files_sort[0]: $!";
			while (<FILE>) {
				chomp;
				my @line = split '_';
				unless (/^[@]/) {
					$split_files_key_read{$split_files_sort[0]} = $line[0];
					print " $split_files_sort[0]: $line[0]\n" if $quiet eq 'no';;
					last;
				}
			}
			if (scalar(keys %split_files_key_read) == $n) {
				print " First reads were recorded successfully.\n" if $quiet eq 'no';
			} else {
				print "Fail to record first reads for $n pieces of FA.\nFatal error. Aborted.\n";
				die "Fail to record first reads for $n pieces of FA: $!";
			}
		} else {
			open FILE, "<", $input_dir.$in_raw or die "cannot open $input_dir${in_raw}: $!";
			while (<FILE>) {
				chomp;
				#my @line = split /\t/;
				if (/^>/) {
					my @line = split '_';
					$split_files_key_read{$in_raw} = $line[0];
					last;
				}
			}
			close FILE;
		}
	}