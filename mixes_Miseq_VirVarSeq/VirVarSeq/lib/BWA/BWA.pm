package BWA;
 
use strict;
use Log;

###### Default BWA parameters

###### Functions

sub alignRead {
	my ($sample, $read, $reference, $nr, $rundir, $params) = @_;
	if(not -e $reference.".bwt"){
	        Log::screen($0, "bwa start indexing reference $reference");
		system("bwa index $reference");
	}
	Log::screen($0, "Start aligning read");
	my $sai = $rundir."/".$sample.".".$nr.".sai";
	system("bwa aln $params -f $sai $reference $read");
	Log::screen($0, "End aligning read");
	return $sai;
}	

sub mapSingle{
	my ($sample, $read1, $read2, $reference, $rundir, $samfile, $params) = @_;
	# Align read 1
	my $sai_1 = alignRead($sample, $read1, $reference, 1, $rundir, $params);
	Log::screen($0, "Start generating alignment for read 1");
	system("bwa samse $reference $sai_1 $read1 > $rundir/$samfile ");
	Log::screen($0, "End generating alignment for read 1");
	
	# If exists, align read2
	if($read2 ne "") {
		my $sai_2 = "";
		$sai_2 = alignRead($sample, $read2, $reference, 2, $rundir, $params);
		Log::screen($0, "Start generating alignment for read 2");
		system("bwa samse $reference $sai_2 $read2 | grep -v -e \"^@\"  >> $rundir/$samfile ");
		Log::screen($0, "End generating alignment for read 2");
	}
}

sub mapPaired{
	my ($sample, $read1, $read2, $reference, $rundir, $samfile, $params) = @_;
	# Align read 1
	my $sai_1 = alignRead($sample, $read1, $reference, 1, $rundir, $params);
	# Align read 2
	my $sai_2 = alignRead($sample, $read2, $reference, 2, $rundir, $params);
	Log::screen($0, "Start generating paired alignment");
	system("bwa sampe -P $reference $sai_1 $sai_2 $read1 $read2 -f $rundir/$samfile");
	Log::screen($0, "End generating paired alignment");
}

1;
