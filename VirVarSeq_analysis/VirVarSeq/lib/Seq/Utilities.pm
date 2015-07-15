package Seq::Utilities; 

use strict;
use warnings;

use Files;
use Log;
use Seq::Constants qw(getIUPAC getAminoAcid getAA);
use Statistics::Basic qw(:all nofill);

our @directions= qw(forward reverse);

my %genetic_code = (
   'TCA' => 'S',    # Serine
   'TCC' => 'S',    # Serine
   'TCG' => 'S',    # Serine
   'TCT' => 'S',    # Serine
   'TTC' => 'F',    # Phenylalanine
   'TTT' => 'F',    # Phenylalanine
   'TTA' => 'L',    # Leucine
   'TTG' => 'L',    # Leucine
   'TAC' => 'Y',    # Tyrosine
   'TAT' => 'Y',    # Tyrosine
   'TAA' => '.',    # Stop
   'TAG' => '.',    # Stop
   'TGC' => 'C',    # Cysteine
   'TGT' => 'C',    # Cysteine
   'TGA' => '.',    # Stop
   'TGG' => 'W',    # Tryptophan
   'CTA' => 'L',    # Leucine
   'CTC' => 'L',    # Leucine
   'CTG' => 'L',    # Leucine
   'CTT' => 'L',    # Leucine
   'CCA' => 'P',    # Proline
   'CCC' => 'P',    # Proline
   'CCG' => 'P',    # Proline
   'CCT' => 'P',    # Proline
   'CAC' => 'H',    # Histidine
   'CAT' => 'H',    # Histidine
   'CAA' => 'Q',    # Glutamine
   'CAG' => 'Q',    # Glutamine
   'CGA' => 'R',    # Arginine
   'CGC' => 'R',    # Arginine
   'CGG' => 'R',    # Arginine
   'CGT' => 'R',    # Arginine
   'ATA' => 'I',    # Isoleucine
   'ATC' => 'I',    # Isoleucine
   'ATT' => 'I',    # Isoleucine
   'ATG' => 'M',    # Methionine
   'ACA' => 'T',    # Threonine
   'ACC' => 'T',    # Threonine
   'ACG' => 'T',    # Threonine
   'ACT' => 'T',    # Threonine
   'AAC' => 'N',    # Asparagine
   'AAT' => 'N',    # Asparagine
   'AAA' => 'K',    # Lysine
   'AAG' => 'K',    # Lysine
   'AGC' => 'S',    # Serine
   'AGT' => 'S',    # Serine
   'AGA' => 'R',    # Arginine
   'AGG' => 'R',    # Arginine
   'GTA' => 'V',    # Valine
   'GTC' => 'V',    # Valine
   'GTG' => 'V',    # Valine
   'GTT' => 'V',    # Valine
   'GCA' => 'A',    # Alanine
   'GCC' => 'A',    # Alanine
   'GCG' => 'A',    # Alanine
   'GCT' => 'A',    # Alanine
   'GAC' => 'D',    # Aspartic Acid
   'GAT' => 'D',    # Aspartic Acid
   'GAA' => 'E',    # Glutamic Acid
   'GAG' => 'E',    # Glutamic Acid
   'GGA' => 'G',    # Glycine
   'GGC' => 'G',    # Glycine
   'GGG' => 'G',    # Glycine
   'GGT' => 'G',    # Glycine
);


sub readReference{
   my ($fasta_file) = @_;
   my @lines = Files::openFile($fasta_file);
   my $reference = "";
   foreach my $line(@lines){
      if(not($line =~ />/)){
         $reference .= $line;
      }
   }
   return $reference;
}

sub translateCodon {
   my($codon) = @_;
   my $newcodon = "";
   my @pos1= @{getIUPAC(substr($codon,0,1))};
   my @pos2= @{getIUPAC(substr($codon,1,1))};
   my @pos3= @{getIUPAC(substr($codon,2,1))};

   #my @pos1= @{$IUPAC_nucleotides{substr($codon,0,1)}};
   #my @pos2= @{$IUPAC_nucleotides{substr($codon,1,1)}};
   #my @pos3= @{$IUPAC_nucleotides{substr($codon,2,1)}};
   for(my $i=0; $i<scalar(@pos1);$i++){
      for(my $j=0; $j<scalar(@pos2);$j++){
         for(my $k=0; $k<scalar(@pos3);$k++){
            $newcodon .= $pos1[$i].$pos2[$j].$pos3[$k].",";
         }
      }
   }
   $newcodon =~ s/,$//;
   return($newcodon);
}

sub codon2Amino {
   # input can be string of codons like CCC,CCT
   # output can be a string of amino acids like P,L
   my($codons) = @_;
   my @codons  = split(/,/, $codons);
   my %amino;
   my $amino_acid;
   foreach my $codon (@codons) {
      #$amino_acid = $AMINO_ACID{$codon};
	  $amino_acid = Seq::Constants::getAA($codon);
      $amino{$amino_acid} = $amino_acid;
   } 
   $amino_acid = join(",", keys %amino);
   return($amino_acid);
}

#
#   Translate a DNA 3-character codon to an amino acid
#
sub codon2aa {
	my ($codon) = @_;
	$codon = uc $codon;
	if (exists $genetic_code{$codon}) {
		return ($genetic_code{$codon});
	} elsif ($codon eq '---'){
		return '';
	} else {
		return 'X';
	}
}
		
#
#   Calculate CodonTable and print to STDOUT
#
sub sam2CodonTable {
	my $sam_file 	= shift;
	my $ref_file 	= shift;
	my $do_trimming = shift;
	my $region_start= shift;
	my $region_len	= shift;
	my $min_qual    = shift;
	my $sample_name = shift;
	my $stats_file  = shift;
	my $codon_file  = shift;
	
	my $ref = readReference($ref_file);

	my %codon_table = ();
	my %codon_table_mins = ();
	my %denom	= ();
	my %codons      = ();

	my $line;
	my $ref_start	    = 1;
	my $indel_reads	    = 0;
	my $min_qual_count  = 0;
	my $count_x	    = 0;
	my $fwd_freq	    = 0;
	my $rev_freq	    = 0;
	my $low_freq	    = 0;
	my $low_freq_dir    = "n/a";
	my %min_values	    = ();
	my $nor  	    = 0; # number or reads
	my $noc  	    = 0; # number or codons
        my $qual;	
	
	open(IN,$sam_file) or die "Could not open $sam_file\n";
	while ($line = <IN>){  
		next if ($line =~ /^@/);
		$nor++;	
		my @input = split /\t/, $line;
		my $seq = $input[0];
		my $direction = $input[1];
	  
		if ($direction & 16){
			$direction = 'reverse';
		} else {
			$direction = 'forward';
		}
	  
		my $pos 		= $input[3]; # start position towards the reference
		next if $pos == 0;  
		my $match 		= $input[5]; # cigar
		my $read 		= $input[9];
		my @nt_arr = split //, $read;
		$qual 		    = $input[10];	
		my @qual_arr    = split //,$qual;
		my $trim_start  = 0;
		my $trim_end    = 0;
		my $start_pos = $ref_start + $pos - 1;

		if ($match =~ m/[ID]/) {
			$indel_reads++;
			next;
		}
	
		if (($match =~ m/^([0-9]+)M/) && ($match =~ m/([0-9]+)M$/) && ($match =~ m/S/)){
			print "SOFT clipping in middle of read not supported in $match !\n";
			exit(1);
		} 
	  
		if ($match=~ m/^([0-9]+)S/){
			$trim_start= $1;
			if ($do_trimming == 0){		# sam file contains start position after trimming.
				$start_pos = $start_pos-$trim_start;
			} elsif ($do_trimming == 1){
				for (my $i=1; $i <= $trim_start; ++$i){
					shift @nt_arr;
					shift @qual_arr;
				}
			}
		}
		
		if (($do_trimming == 1) && ($match=~ m/([0-9]+)S$/)){
			$trim_end= $1;
			for (my $i=1; $i <= $trim_end; ++$i){
				pop @nt_arr;
				pop @qual_arr;
			}
		}

		my %nt_table   = ();
		my %qual_table = ();
		my $prev_pos   = $start_pos - 1;
		while ($match =~ m/([0-9]+)([MS])/g){
			my $len = $1;
			my $type = $2;	
			if (($type eq 'S') && ($do_trimming == 1)) {
				# do nothing
			} else {	
				my $pos= $prev_pos;			
				for (my $i=1; $i <= $len; ++$i){				
					$pos= $prev_pos + $i;
					my $nt= "-";
					$qual = "~";
					$nt = shift @nt_arr;
					$qual = shift @qual_arr;
					$nt_table{$pos} = $nt;
					$qual_table{$pos} = $qual;
				}
				$prev_pos= $pos;					
			}				
		}
		my $codon_len= 0;
		my $codon= "";
		my $orf= -1;
		my $num_codons= 0;  

		# $pos = foreach nucleotide pos
		# orf = open reading frame, number of nucleotides missing at beginning of read
		# num_codons = codon number of reference (not number of codons!)
		# this foreach loop only specifies for the first codon of the entire read the number of missing nucleotides. These missing nucleotides are filed with "-" (see for loop below)
		foreach my $pos (sort {$a <=> $b} keys %nt_table){
			if ($pos >= $region_start){
				$orf= ($pos - $region_start)%3;			
				$num_codons= int(($pos - $region_start)/3);
				last;		
			}			
		}
	
		for (my $i=1; $i <= $orf; ++$i){
			$codon .= "-";
			$codon_len += 1;
		}   
   	
		my $last_pos  = 0;
		my $score     = 0;
		$qual	  = 0;   
		my $min_value = 0;
		my @scores;	
		my @codon_means;
		my @codon_mins;

		
		foreach my $pos (sort {$a <=> $b} keys %nt_table){
			$last_pos = $pos;		
			if ($pos >= $region_start){
				my $nt = $nt_table{$pos};				
				$codon .= $nt;
				$codon_len += 1;
				$qual = $qual_table{$pos};			
				$score =  ord($qual) - 33;
				push(@scores, $score);
				if ($codon_len == 3){
					$num_codons += 1;
					if ($num_codons > $region_len){
						last;
					}
					if ($orf > 0){
						$orf= 0;
					} else {
						my @sorted_scores = sort { $a <=> $b } @scores;
						$min_value = $sorted_scores[0];
						
						if (exists($min_values{$min_value})) {
							$min_values{$min_value}++;
						}else{
							$min_values{$min_value}=1;
						}						
						
						$noc++;
						if ($min_value < $min_qual) {   # discard all codons of which the mininum quality is lower then the min_qual value
							$min_qual_count++;
						} else {			
							$codon_table{$num_codons}{$codon}{$direction} += 1;
							$denom{$num_codons}{$direction} += 1;
							$codons{$codon}= 1;
							push( @{$codon_table_mins{$num_codons.$codon.$direction}}, $min_value );				
						}
					}
					$codon       = "";
					$codon_len   = 0;
					@scores      = ();
					@codon_means = ();
					@codon_mins  = ();
				}
			}
		}
   	
		if (($codon_len > 0) && ($num_codons <= $region_len)){
			my $keep = 1;
			if ($codon_len == 1){
				$codon .= "--";
				$keep= 0;
			}
			if ($codon_len == 2){
				$codon .= "-";
				$keep= 0;
			}
			$num_codons += 1;
			if ($keep == 1){
				$codon_table{$num_codons}{$codon}{$direction} += 1;
				$denom{$num_codons}{$direction} += 1;
				$codons{$codon}= 1;
			}
		}
	}
	close(IN);
   
	open(CODON,">$codon_file") or die "Could not open [$codon_file] for writing codon table\n";   
   
	print CODON "SAMPLE\tPOSITION\tREF_CODON\tCODON\tREF_AA\tAA\tFWD_CNT\tFWD_DENOM\tREV_CNT\tREV_DENOM\tFWD_MEAN_MIN_QUAL\tREV_MEAN_MIN_QUAL\tFWD_FREQ\tREV_FREQ\tLOW_FREQ\tLOW_FREQ_DIRECTION\tFWD_STDDEV_MIN_QUAL\tREV_STDDEV_MIN_QUAL\tCNT\tDENOM\tFREQ\n";   
   
	# This $pos = $codon_number, not nucleotide pos	
	foreach my $pos (sort {$a <=> $b} keys %codon_table){
		my $ins_flag = 0;
		if ($pos =~ m/\./){
			$ins_flag = 1;
		}
		my $ref_codon = "---";
		if ($ins_flag == 0){
			$ref_codon = substr($ref,($pos-1)*3+$region_start-$ref_start,3);
		}
		my %cnt2 = ();
		if ($ins_flag == 1){
			foreach my $dir (@directions){
				foreach my $ins_codon (sort keys %{$codon_table{$pos}}){
					if ($ins_codon ne '---'){
						$cnt2{$dir} += $codon_table{$pos}{$ins_codon}{$dir};
					}
				}
			}
		}
   
		foreach my $codon (sort keys %codons){   				
			my %cnts = ();
			my %tots = ();
			my $rev_mean_min_qual   = 0;
			my $fwd_mean_min_qual   = 0;	
			my $rev_stddev_min_qual = 0;	
			my $fwd_stddev_min_qual = 0;	
			my $coverage            = 0;	
			my $cnts_codon          = 0;
			my $frequency           = 0;
   						
			my $prefix = "";
			if ($pos < 10){
				$prefix .= "0";
			}
			if ($pos < 100){
				$prefix .= "0";
			}													
						
			foreach my $dir (@directions){
				my $cnt = 0;
				my $tot = 0;
				my @median_min_qual = 0;
				if (exists $denom{$pos}{$dir}){
					$tot= $denom{$pos}{$dir};
				}
				if ($ins_flag == 1){
					my $pos2= $pos;
					$pos2=~ s/\.[0-9]+$//;
					$tot= $denom{$pos2}{$dir};
				}
				if (exists $codon_table{$pos}{$codon}{$dir}){
					$cnt= $codon_table{$pos}{$codon}{$dir};   
					if ($dir eq "reverse") {
						$rev_mean_min_qual = mean(@{$codon_table_mins{$pos.$codon.$dir}});
						$rev_stddev_min_qual = stddev(@{$codon_table_mins{$pos.$codon.$dir}});					
					}
					if ($dir eq "forward") {
						$fwd_mean_min_qual = mean(@{$codon_table_mins{$pos.$codon.$dir}});
						$fwd_stddev_min_qual = stddev(@{$codon_table_mins{$pos.$codon.$dir}});										
					}   				
				}
				if (($ins_flag == 1) && ($codon eq '---')){
					$cnt= $tot - $cnt2{$dir};
				}
				$cnts{$dir}= $cnt;
				$tots{$dir}= $tot;
			}
			if (($cnts{'forward'} > 0) || ($cnts{'reverse'} > 0)){
				print CODON "$sample_name\t$prefix$pos\t$ref_codon";
				my $ref_aa = codon2aa($ref_codon);			
				my $aa = codon2aa($codon);
				if ($aa eq "X"){
					$count_x = $count_x + $cnts{'forward'} + $cnts{'reverse'};
				}
				if (($ins_flag == 1) && ($codon ne '---')){
					$codon = "+$codon";
					$aa = "+$aa";
				}
				print CODON "\t$codon\t$ref_aa\t$aa";	
				if ($tots{'forward'} == 0) {
					$fwd_freq = 0;
				} else {				
					$fwd_freq = $cnts{'forward'}/$tots{'forward'};
				}				
				if ($tots{'reverse'} == 0) {
					$rev_freq = 0;
				} else {							
					$rev_freq = $cnts{'reverse'}/$tots{'reverse'};			
				}				
				if ($fwd_freq < $rev_freq) {
					$low_freq = $fwd_freq;
					$low_freq_dir = "forward";
				} else {
					$low_freq = $rev_freq;
					$low_freq_dir = "reverse";				
				}
				$coverage = $tots{'reverse'} + $tots{'forward'}; # counts pos
				$cnts_codon = $cnts{'reverse'} + $cnts{'forward'}; 
				$frequency = $cnts_codon/$coverage;
				print CODON "\t$cnts{'forward'}\t$tots{'forward'}\t$cnts{'reverse'}\t$tots{'reverse'}\t$fwd_mean_min_qual\t$rev_mean_min_qual\t$fwd_freq\t$rev_freq\t$low_freq\t$low_freq_dir\t$fwd_stddev_min_qual\t$rev_stddev_min_qual\t$cnts_codon\t$coverage\t$frequency\n";
			}
   
			if ($ins_flag == 0){
				my $coverage= $cnts{'forward'} + $cnts{'reverse'};				
			}
		}
	}
    close(CODON);
	
	
	open(STAT,">$stats_file") or die "Could not open [$stats_file] for writing statistics\n";   	
	print STAT "number of reads[$nor]\n";
	print STAT "number of codons[$noc]\n";
	print STAT "number of indels [$indel_reads]\n";
	print STAT "indels [" . ($indel_reads/$nor)*100 . "%]\n";
	print STAT "min_qual_count [$min_qual_count] [" . ($min_qual_count/$noc)*100 ."%]\n";
	print STAT "count_x [$count_x]\n";	
	foreach my $key (sort {$a <=> $b} keys %min_values) {
		my $value = $min_values{$key};
		print STAT "value [$value] key [$key]\n";
	}
	close(STAT);
}
	

#
#   Calculate IndelTable
#
sub sam2IndelTable {
	my $sam_file 	= shift;
	my $ref_file 	= shift;
	my $do_trimming = shift;
	my $region_start= shift;
	my $region_len	= shift;
	my $min_qual    = shift;
	my $sample_name = shift;
	my $indel_file  = shift;	
	
	my $ref = readReference($ref_file);
	
	my %codon_table 		= ();
	my %codon_qual_table 	= ();
	my %codon_table_mins    = ();
	my %denom				= ();
	my $header				= 0;
	my $line;
	my $ref_start  			= 1;
	my $DEBUG				= 1;
	my $indel_reads			= 0;
	my $min_qual_count		= 0;
	my $count_x				= 0;
	my $fwd_freq			= 0;
	my $rev_freq			= 0;
	my $low_freq			= 0;
	my $low_freq_dir		= "n/a";
	my @min_values;
	my %nt_table_total      = ();
	my %nt_table_indel_total= ();

	my $nor 				= 0 ; # number or reads
	my $noc 				= 0; # number or codons	

	open(IN,$sam_file) or die "Could not open $sam_file\n";
	while ($line = <IN>){
		next if ($line =~ /^@/);
		$nor++;	
		my @input 			= split /\t/, $line;
		my $seq				= $input[0];
		my $direction		= $input[1];
		if ($direction & 16){
			$direction = 'reverse';
		}else{
			$direction = 'forward';
		}
		my $pos				= $input[3];   
		my $match 			= $input[5]; #cigar
		my $read			= $input[9];
		my @nt_arr			= split //, $read;
		my @nt_arr_indel	= split //, $read;	
		my $qual			= $input[10];	
		my @qual_arr		= split //,$qual;
		my $trim_start		= 0;
		my $trim_end		= 0;
		my $start_pos		= $ref_start + $pos - 1; 

		if ($match =~ m/[ID]/) {
			if ($match=~ m/^([0-9]+)I/){
				print "INS at start of read not supported in $match !\n";
				exit(1);
			}
			if ($match=~ m/^([0-9]+)D/){
				print "DEL at start of read not supported in $match !\n";
				exit(1);
			}
			if ($match=~ m/([0-9]+)I$/){
				print "INS at end of read not supported in $match !\n";
				exit(1);
			}
			if ($match=~ m/([0-9]+)D$/){
				print "DEL at end of read not supported in $match !\n";
				exit(1);
			}
			$indel_reads++;
		}	
	
		if (($match=~ m/^([0-9]+)M/) && ($match=~ m/([0-9]+)M$/) && ($match=~ m/S/)){
			print "SOFT clipping in middle of read not supported in $match !\n";
			exit(1);
		}
		if ($match=~ m/^([0-9]+)S/){
			$trim_start= $1;
			if ($do_trimming == 0){		#sam file contains start position after trimming.
				$start_pos= $start_pos-$trim_start;
			}
			elsif ($do_trimming == 1){
				for (my $i=1; $i <= $trim_start; ++$i){
					shift @nt_arr;
					shift @qual_arr;
					shift @nt_arr_indel;		
				}
			}
		}
		if (($do_trimming == 1) && ($match=~ m/([0-9]+)S$/)){
			$trim_end= $1;
			for (my $i=1; $i <= $trim_end; ++$i){
				pop @nt_arr;
				pop @qual_arr;
				pop @nt_arr_indel;
			}
		}

		my %nt_table= ();
		my %nt_table_indel= ();
		my $prev_pos= $start_pos - 1;			
		while ($match =~ m/([0-9]+)([MIDS])/g){
			my $len= $1;
			my $type= $2;
			if (($type eq 'M') || (($type eq 'S') && ($do_trimming == 0)) || ($type eq 'D')){
				my $pos= $prev_pos;
				for (my $i=1; $i <= $len; ++$i){
					$pos= $prev_pos + $i;
					my $nt= "-";
					if ($type eq 'D'){
						$nt_table_indel{$pos} = $nt;
						$nt_table{$pos} = $nt;
					}
					else{
						$nt= shift @nt_arr_indel;
						$nt_table{$pos} = $nt;					
					}
				}
				$prev_pos= $pos;
			}
			elsif ($type eq 'I'){
				my $pos= $prev_pos;
				my $ins= "+";
				for (my $i=1; $i <= $len; ++$i){
					my $nt= shift @nt_arr_indel;
					$ins .= "$nt";
				}
				$nt_table_indel{$pos}= $ins;			
			}
			else{
				#do nothing.
			}
		}
	
		while ( my($nt_table_indel_pos, $nucl) = each(%nt_table_indel) ) {
			$nt_table_indel_total{$nt_table_indel_pos.";".$nucl.";".$direction} += 1;
		}

		while ( my($nt_table_pos, $nt_table_value) = each(%nt_table) ) {
			$nt_table_total{$nt_table_pos.";".$direction} += 1;
		}
	}
	close(IN);

	#
	# Print INDELS
	#
	my $max_indel_freq = 0;
	open(INDEL,"> $indel_file") or die "ERROR: cannot open [$indel_file]\n";
	print INDEL "POS\tINDEL\tDIRECTION\tCOUNT\tTOTAL COUNT\tFREQ\n";
	foreach my $nt_table_indel_total_pos (sort {$a cmp $b} keys %nt_table_indel_total) {
		my $indel_count     = $nt_table_indel_total{$nt_table_indel_total_pos};
		my @pos_values      = split(/;/, $nt_table_indel_total_pos);
		my $indel_pos       = $pos_values[0];
		my $indel_indel     = $pos_values[1];
		my $indel_direction = $pos_values[2];     
		my $total_count     = $nt_table_total{$indel_pos.";".$indel_direction}; # denominator prev pos (forward and reverse)
		my $freq = $indel_count/$total_count;   
		if ($freq > $max_indel_freq) {
			$max_indel_freq = $freq;
		}
		print INDEL "$indel_pos\t$indel_indel\t$indel_direction\t$indel_count\t$total_count\t$freq\n";
	}
	close(INDEL);
		
}	


	
1;
