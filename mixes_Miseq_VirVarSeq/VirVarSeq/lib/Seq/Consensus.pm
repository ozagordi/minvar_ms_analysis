package Seq::Consensus;

use strict;
use warnings;

my $READLENGTH = "";


sub getConsensus{
	my ($mapping, $reference, $consensus_start, $consensus_end, $ROIstart, $ROIend, $do_trimming, $correct_consensus)= @_;
	
	# first read in full sam file
	my $nucl_table_ref = readSamFile($mapping, $reference, $do_trimming, $consensus_start, $consensus_end, $ROIstart, $ROIend);

	# then write out consensus
	my $consensus = buildConsensus($nucl_table_ref, $reference, $consensus_start, $consensus_end, $correct_consensus, $ROIstart, $ROIend);
	return $consensus;
}

sub readSamFile{
	my ($mapping, $reference, $do_trimming, $consensus_start, $consensus_end, $ROIstart, $ROIend)= @_;

	my $ref_start = 1;	
	my %nucl_table =();
	my $read_length;
	
   	# Loop over all reads in sam file
	open(IN,$mapping) or die "Could not open $mapping\n";
	while (<IN>){
		$_ =~ s/\r|\n//sg;
		if(not($_ =~ /^\@/)){
			my @samstring= split("\t", $_);
			if (scalar @samstring < 10){
				last;
			}
			my $seq 	= $samstring[0];
			my $samflag = $samstring[1];
			my $pos		= $samstring[3];
			my $match	= $samstring[5];
			my $read	= $samstring[9];
			
			my $direction = "";
			if ($samflag & 16){
				$direction = 'reverse';
			}else{
				$direction = 'forward';
			}
			my $trim_start	= 0;
			my $trim_end	= 0;
			my $consensus_start = $ref_start+$pos-1;			 
			
			
			# only change global read length parameter when longer
			$read_length = setReadLength(length($read));
			 
			# ignore all reads with insertions or deletions at begin or end and softclipping in the middle
			if ($match=~ m/^([0-9]+)I/ or $match=~ m/^([0-9]+)D/
				or $match=~ m/([0-9]+)I$/ or $match=~ m/([0-9]+)D$/
				or (($match=~ m/^([0-9]+)M/) &&
				($match=~ m/([0-9]+)M$/) && ($match=~ m/S/))){
                                # do nothing
			}
			else{
				my @nt_arr= split(//, $read);
				if ($match=~ m/^([0-9]+)S/){
					$trim_start= $1;
					if ($do_trimming == 0){      
						#sam file contains start position after trimming.
						$consensus_start= $consensus_start-$trim_start;
					}
					elsif ($do_trimming == 1){
						for (my $i=1; $i <= $trim_start; ++$i){
							shift @nt_arr;
						}
					}
				}
				elsif (($do_trimming == 1) && ($match=~ m/([0-9]+)S$/)){
					$trim_end= $1;
					for (my $i=1; $i <= $trim_end; ++$i){
						pop @nt_arr;
					}
				}

				my $prev_pos= $consensus_start - 1;
				while ($match =~ m/([0-9]+)([MIDS])/g){
					my $len= $1;
					my $type= $2;
					if (($type eq 'M') || (($type eq 'S') && ($do_trimming == 0)) || ($type eq 'D')){
						my $pos= $prev_pos;
						for (my $i=1; $i <= $len; ++$i){
							$pos= $prev_pos + $i;
							my $nt= "-";
							if ($type eq 'D'){
								#do nothing.
							}
							else{
								$nt= shift @nt_arr;
							}
							$nucl_table{$pos}{$nt}{$direction} += 1;
						}
						$prev_pos= $pos;
					}
					elsif ($type eq 'I'){
						my $pos= $prev_pos;
						my $ins= "+";
						for (my $i=1; $i <= $len; ++$i){
							my $nt= shift @nt_arr;
							$ins .= "$nt";
						}
						$nucl_table{"$pos.1"}{$ins}{$direction} += 1;
					}
					###else{
					###	#do nothing.
					###}
				}
			}
		}
	}
	close(IN);
	return \%nucl_table;
}

sub buildConsensus{
	my ($nucl_table_ref, $reference, $consensus_start, $consensus_end, $correct_consensus, $ROIstart, $ROIend)= @_;
	
	my %nucl_table = %$nucl_table_ref;
	my $consensus = "";
	my $consensus_type = "";
	my %inserts =();
	my @directions= ("forward", "reverse");
	my $read_length = getReadLength();
	
	my %nts= ();
	$nts{'A'}= 1;
	$nts{'C'}= 1;
	$nts{'G'}= 1;
	$nts{'T'}= 1;
	$nts{'N'}= 1;
	$nts{'-'}= 1;		
	
	for (my $i=0; $i < length($reference); ++$i){
		my $pos= 1 + $i;
		my $ref_nt= substr($reference,$pos-1,1);
		my $cons_nt= $ref_nt;
		my $max_cnt= 0;
		if (($pos >= $consensus_start) && ($pos <= $consensus_end)){		
			foreach my $nt (sort keys %nts){
				my $cnt= 0;
				foreach my $dir (@directions){
					if (exists $nucl_table{$pos}{$nt}{$dir}){
						$cnt += $nucl_table{$pos}{$nt}{$dir};           					
					}
				}
				if ($cnt > $max_cnt){
					$cons_nt= $nt;
					$max_cnt= $cnt;
				}
			}
		}
		
		if ($cons_nt eq "-") {
			$cons_nt = "N";
			$consensus_type .= "D";
		}else{
			$consensus_type .= "M";
		}
		$consensus .= $cons_nt;			
		
		#
		# Keep track of Inserts
		#	  
		my $coverage = 0;
		my $coverage_ins = 0;
		if ($correct_consensus == 1 && $ROIstart <= $pos && $ROIend >= $pos) {  
			$max_cnt = 0;	
			if (($pos >= $consensus_start) && ($pos <= $consensus_end)){
				foreach my $nt (sort keys %nucl_table){
					foreach my $dir (@directions){						
						if (exists $nucl_table{"$pos.1"}{$nt}{$dir}){							
							$coverage_ins += $nucl_table{"$pos.1"}{$nt}{$dir};  
						}
						if (exists $nucl_table{$pos}{$nt}{$dir}){							
							$coverage += $nucl_table{$pos}{$nt}{$dir};                     
						}
					}
				}
				$max_cnt = $coverage - $coverage_ins;
				
				$cons_nt = "I";
				foreach my $nt (sort keys %nucl_table){
					my $cnt= 0;
					foreach my $dir (@directions){						
						if (exists $nucl_table{"$pos.1"}{$nt}{$dir}){
							$cnt += $nucl_table{"$pos.1"}{$nt}{$dir};  
						}
					}			
					if ($cnt > $max_cnt){
						$cons_nt= $nt;
						$max_cnt= $cnt;
					}			
				}		
				if ($cons_nt ne "I") {
					$inserts{$pos} = $cons_nt;
				}
			}			
		}
	}
	
	while (my($insert_pos, $insert_nt) = each(%inserts)){		
		my $region_to_check  = substr($consensus_type, ($insert_pos - $read_length), ($read_length * 2));
		if ($region_to_check =~ m/D{1}/) {  # only one D exists 		
			#
			# HANDLE DELETE
			#
			# on deleted pos => remove base in consensus
			# ==> check to prevent 1 D otherwise error					
			my $d_pos = index($region_to_check, "D");											
			my $consensus_pos = $d_pos + $insert_pos - $read_length;					
			substr($consensus, $consensus_pos, 1, '');				
			#
			# HANDLE INSERT
			#			
			# on $insert_pos => add base in consensus (no replace) !!!!		
		        my $consensus_before_insert = substr($consensus, 0, $insert_pos-1);
			my $consensus_after_insert  = substr($consensus, $insert_pos-1);
			# strip the + sign
			$insert_nt =~ s/\+//;
			$consensus = $consensus_before_insert . $insert_nt . $consensus_after_insert;			
		}else{
			# D can stay in consensus_type(= N in consensus)
			print "!!! WARNING : 0 or more then 1 DELETION found during consensus correction\n";
			print "!!! region_to_check [$region_to_check], insert_pos [$insert_pos], insert_nt [$insert_nt], consensus_type [$consensus_type]\n";
		}
	}	
	
	return($consensus);	
}

sub getReadLength{
	return $READLENGTH;
}

sub setReadLength{
	my ($len) = @_;
	if ( $READLENGTH eq "" or $len > $READLENGTH) {	$READLENGTH = $len; }
	return $len;
}


1;
