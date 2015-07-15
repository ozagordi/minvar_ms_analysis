#!/usr/bin/env perl

use strict;
use warnings;

use Files;
use Log;
use Seq::Consensus;
use Seq::Utilities;
use Getopt::Long;
use Pod::Usage;

### Set default values or initiate variables #######################################
my ($sample_file, $ref, $project_in, $project_out, $optionsfile)  = "";
my $consensus_start 	= 0;
my $consensus_end 		= 0;
my $ROIstart 			= 0;
my $ROIend 				= 0;
my ($man, $help) 		= 0;
my $trimming			= 1;
my $correct_consensus	= 0;
my $map_vs_ref			= "map_vs_ref";
my $consensus			= "consensus";

### Read options from command line #################################################
GetOptions(
	'ref=s'			=> \$ref, 
	'samplelist=s'		=> \$sample_file ,	
	'outdir=s'		=> \$project_out ,
	'indir:s'		=> \$project_in , 
	'start=i'		=> \$consensus_start,
	'end=i'			=> \$consensus_end, 
	'ROIstart:i'		=> \$ROIstart,
	'ROIend:i'		=> \$ROIend,
	'trimming:i'		=> \$trimming,
	'correct:i'		=> \$correct_consensus,
	'refmapdir:s'		=> \$map_vs_ref,
	'consensusdir:s'	=> \$consensus,
	'help!'			=> \$help ,		
	'man!'			=> \$man		#,
	#'options=s'	=>  \$optionsfile
) or pod2usage(2); 
pod2usage("\n") if $help;
pod2usage(-verbose => 2) if $man;

### 2.0 Check command line arguments ###################################################
Log::screen($0, "2.0 Checking command line arguments ...");
checkParameters();

### 2.1-2.4 Calculate consensus sequences ##################################################

### Set file and directory locations
Log::screen($0, "2.1 Checking input and output locations ...");
#my ($map_vs_ref, $consensus) = 
checkFileLoc();

# Check references
Log::screen($0, "2.2 Load in reference ...");
my $reference = Seq::Utilities::readReference($ref);

# Read in sample files
Log::screen($0, "2.3 Get sample names from file [$sample_file]...");
my @samples = Files::openFile($sample_file);
my $num_of_samples = scalar(@samples);

# Calculate consensus file for each sample
Log::screen($0, "2.4 Calculate consensus ...");
for (my $i=0; $i < $num_of_samples; ++$i){
	my $sample = $samples[$i];
	Log::screen($0, "-> for sample [$sample] ...");
	
	my $sam_file = $project_out."/".$map_vs_ref."/".$sample.".sam";
	my ($consensusseq) = Seq::Consensus::getConsensus($sam_file, $reference, $consensus_start, $consensus_end, $ROIstart, $ROIend, $trimming, $correct_consensus);
	#1, 1, "$map_vs_ref/$sample.sam", $start_pos, $end_pos, $reference);
	my $out_file = $project_out."/".$consensus."/".$sample."_consensus.fa";
	open(CONSENSUS,">$out_file") or die "ERROR: cannot open file [$out_file] for writing consensus: $!\n";
	print CONSENSUS ">".$sample."\n".$consensusseq."\n";
	close(CONSENSUS);
}
Log::screen($0, "2.5 End.");

########################################################################################
### AUX FUNCTIONS
########################################################################################

sub checkParameters{
	#if($optionsfile ne ""){ 
	#	parseOptionsFile();
	#}
	#else {
	### Check for required arguments 
	if($sample_file eq ""){
	   pod2usage("sample file (sample=)" );
	   exit(1);
	}
	if($project_out eq ""){
	   pod2usage("output folder for project (outdir=)" );
	   exit(1);
	}
	if($consensus_start == 0){
	   pod2usage("start position (start=)" );
	   exit(1);
	}
	if($consensus_end == 0){
	   pod2usage("end position (end=)" );
	   exit(1);
	}
	if($project_out eq ""){
	   pod2usage("reads folder for project (indir=)" );
	   exit(1);
	}
	if($ref eq ""){
	   pod2usage("reference (ref=)" );
	   exit(1);
	}
	#}
	# Additional check for region of interest parameters
        if ( $ROIstart eq 0 ) {
           $ROIstart = $consensus_start;
	}	
        if ( $ROIend eq 0 ) {
           $ROIend = $consensus_end;
	}	
}

#######################
#sub ParseOptionsFile{
#
#}

#######################
sub checkFileLoc{
	#	Set location of input and output files
	my $map_vs_refdir	  = $project_out."/".$map_vs_ref;
	if (not -d $map_vs_refdir) {
		die "ERROR: Directory with reference mappings does not exist.";
	}
	# Create folder for storing consensus fasta files if not present	
	my $consensus_dir = $project_out."/".$consensus;
	Log::screen($0, "|--- Create local [$consensus_dir] dir ...");
	unless(-d $consensus_dir){
		mkdir $consensus_dir or die "ERROR: Cannot create directory for consensus files." ;
	}
}
	
#######################
sub printUsage {
   print "\n\n";
   print "\tusage : consensus.pl <source> <target> <sample_file>\n";
   print "\t-----\n";
   print "\twith\n";
   print "\t\t<sample_file> a file with the name of the fastq samples (from SampleSheet.xls) to calculate the consensus of \n";
   print "\t\t<ref> name of the reference on the CLCBio Server\n";
   print "\t\t   possible values are\n\n";
   #while ( my ($key, $value) = each(%references) ) {
   #   print "\t\t\t$key \t => $value\n";
   #}
   print "\n\n";
   print "\t\t<start_pos> start position\n";
   print "\t\t<end_pos> end position\n";
   print "\n\n";
   print "\texample\n";
   print "\t\tconsensus.pl samples.txt JCV_AB038249 2012 3879\n";
   print "\n\n";
   print "\n\n";
}

#######################
exit(0);


__END__
=head1 NAME

consensus.pl 
The consensus at each position of the reference genome will be determined.

=head1 SYNOPSIS

perl consensus.pl --samplelist <file> --ref <ref> --outdir <outdir> [options]

Options:

--help brief help message

--man full documentation

=head1 OPTIONS

=over 8

=item --help

Program:	consensus.pl 

Version:	2013-10-01

Contact:	Joke Reumers (jreumers@its.jnj.com); Yves Wetzels (ywetzel@its.jnj.com)

Usage:   	perl consensus.pl [options] 

Options:
	'ref=s'			=> \$ref, 
	'samplelist=s'		=> \$sample_file ,	
	'outdir=s'		=> \$project_out ,
	'indir:s'		=> \$project_in , 
	'start=i'		=> \$consensus_start,
	'end=i'			=> \$consensus_end, 
	'ROIstart:i'		=> \$ROIstart,
	'ROIend:i'		=> \$ROIend,
	'trimming:i'		=> \$trimming,
	'correct:i'		=> \$correct_consensus,
	'refmapdir:s'		=> \$map_vs_ref,
	'consensusdir:s'	=> \$consensus,
	'help!'			=> \$help ,		
	'man!'			=> \$man	


	--samplelist		<string>	[ required ]
	--ref			<string>	[ required ]
	--start			<integer>	[ required ]
	--end			<integer>	[ required ]
	--indir       		<string>	[ required ]
	--outdir		<string>	[ required ]	
	
	--ROIstart          	<integer>	[ optional] [ default = amplicon start ]
	--ROIend            	<integer>	[ optional] [ default = amplicon end ]
	--trimming          	<integer>	[ optional] [ default: 1 = trimming on; 0 = trimming off ]	
	--correct		<integer>	[ optional] [ default: 1 = correct consensus on; 0 = off ]	
	--refmapdir		<string>	[ optional] [ default: map_vs_ref ]
	--consensusdir		<string>	[ optional] [ default: consensus ]

For details on these options run 
	perl consensus.pl --man
	
=item --man

	-s/--samplelist	<samplelist>    Text file containing the sample names as used in the fastq filenames.	
	-r/--ref <reference file>       Full path of the reference file in fasta format. BWA generated indexes should be in the same folder.
	-i/--indir                      Full path of the directory containing the demultiplexed fastq files (in gz format)
	-o/--outdir                     Full path of the target directory for this project


=back

=head1 DESCRIPTION

The consensus at each position of the reference genome will be determined.

=cut
